// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides jst traversal algorithm for algorithms with state.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <libcontrib/closure_adaptor.hpp>
#include <libcontrib/type_traits.hpp>

#include <libjst/container/concept_jst.hpp>
#include <libjst/traversal/concept_execution.hpp>
#include <libjst/traversal/lazy_tree.hpp>
#include <libjst/concept.hpp>

namespace libjst
{
    // Ok now we define the first algorithm:
    // it is a searcher factory creating new searchers

    namespace _searcher_factory
    {
        template <typename jst_node_t, typename sequence_iterator_t>
        struct result
        {
            jst_node_t const & node;
            sequence_iterator_t position;
        };

        template <typename jst_t, typename pattern_t, typename publisher_t>
        struct  stateless_pattern_operation
        {
            jst_t _jst;
            pattern_t _pattern;
            publisher_t _publisher;

        private:

            template <observable_stack subscriber_t>
            constexpr auto make_tree([[maybe_unused]] subscriber_t &subscriber)
            {
                if constexpr (libjst::is_resumable_v<std::remove_cvref_t<pattern_t>>)
                    return lazy_tree{_jst, libjst::window_size(_pattern), subscriber};
                else
                    return lazy_tree{_jst, libjst::window_size(_pattern)};
            }

            // template <typename node_t>
            // constexpr auto label(node_t const & node) const
            // {
            //     auto && seq = node.sequence();
            //     auto it = seq.begin();
            //     size_t head_position{};

            //     if constexpr (!libjst::is_resumable_v<std::remove_cvref_t<pattern_t>>)
            //     {
            //         head_position =  node.first_position() - std::min(libjst::window_size(_pattern),
            //                                                           node.first_position());
            //         assert(seq.size() >= head_position);
            //         assert(node.next_position() >= head_position);
            //         it += head_position;
            //     }
            //     size_t subrange_size = std::min(std::min(node.next_position(), node.last_position()), seq.size()) -
            //                            head_position;
            //     return std::ranges::subrange{it, it + subrange_size, subrange_size}; // borrowed range
            // }

            constexpr friend void tag_invoke(std::tag_t<libjst::start>, stateless_pattern_operation & me) noexcept
            {
                using algorithm_stack_t =
                        std::stack<std::remove_cvref_t<pattern_t>, std::vector<std::remove_cvref_t<pattern_t>>>;

                try {
                    algorithm_stack_t algorithm_stack{};
                    algorithm_stack.push((pattern_t &&)me._pattern);

                    auto jst_tree = me.make_tree(algorithm_stack);

                    for (auto && node : jst_tree) {
                        std::invoke(algorithm_stack.top(), node.sequence(), [&] (auto && it) {
                            libjst::set_next(me._publisher, result{node, it});
                        });
                    }

                    std::cout << "Number of prunes: " << jst_tree._prune_count
                              << " branches: " << jst_tree._branch_count
                              << " total: " << (jst_tree._branch_count + jst_tree._prune_count) << '\n';

                    // finished: so we call finished.
                    libjst::set_value((publisher_t &&)(me._publisher));
                } catch (...) {
                    libjst::set_error((publisher_t &&)(me._publisher), std::current_exception());
                }
            }
        };

        template <typename text_t, typename pattern_t>
        class searcher
        {
        private:
            text_t _text;
            pattern_t _pattern;

        public:

            // constructor
            template <typename _text_t, typename _pattern_t>
            constexpr explicit searcher(_text_t &&text, _pattern_t &&pattern)
                noexcept(std::is_nothrow_constructible_v<text_t, _text_t> &&
                         std::is_nothrow_constructible_v<pattern_t, _pattern_t>) :
                _text{(_text_t &&)text},
                _pattern{(_pattern_t &&)pattern}
            {}
        private:

            template <typename searcher_t, typename publisher_t>
            using operation_t = stateless_pattern_operation<jst::contrib::member_type_t<searcher_t, std::remove_cvref_t<text_t>>,
                                                            jst::contrib::member_type_t<searcher_t, std::remove_cvref_t<pattern_t>>,
                                                            publisher_t>;

            template <typename searcher_t, typename publisher_t>
                requires (std::same_as<std::remove_cvref_t<searcher_t>, searcher>)
                        //   sender_to<std::remove_cvref_t<searcher_t>, publisher_t>)
            constexpr friend auto tag_invoke(std::tag_t<libjst::connect>, searcher_t && me, publisher_t && publisher)
                -> operation_t<searcher_t, publisher_t>
            {
                using fwd_jst_t = jst::contrib::member_type_t<searcher_t,     std::remove_cvref_t<text_t>>;
                using fwd_pattern_t = jst::contrib::member_type_t<searcher_t, std::remove_cvref_t<pattern_t>>;
                return operation_t<searcher_t, publisher_t>{(fwd_jst_t &&)me._text,
                                                            (fwd_pattern_t &&)me._pattern,
                                                            (publisher_t &&)publisher};
            }
        };

        inline constexpr struct _cpo
        {
            template <traversable_journaled_sequence_tree text_t, typename pattern_t>
            constexpr auto operator()(text_t && text, pattern_t && pattern) const
                noexcept(std::is_nothrow_constructible_v<_searcher_factory::searcher<text_t, pattern_t>, text_t, pattern_t>)
                -> _searcher_factory::searcher<text_t, pattern_t>
            {
                return _searcher_factory::searcher<text_t, pattern_t>{(text_t &&)text, (pattern_t &&)pattern};
            }

            template <typename pattern_t>
            constexpr auto operator()(pattern_t && pattern) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, pattern_t>)
                -> jst::contrib::closure_result_t<_cpo, pattern_t>
            { // we need to store the type that needs to be called later!
                return jst::contrib::make_closure(_cpo{}, (pattern_t &&)pattern);
            }
        } jst_searcher{};
    } // namespace _searcher_factory

    using _searcher_factory::jst_searcher;
} // namespace libjst
