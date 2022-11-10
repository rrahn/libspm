// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides jst traversal property for propertys with state.
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
    // Ok now we define the first property:
    // it is a searcher factory creating new searchers

    namespace _searcher_on_lazy_tree
    {
        template <typename jst_node_t, typename sequence_iterator_t>
        struct result
        {
            jst_node_t const & node;
            sequence_iterator_t position;
        };

        template <typename jst_t, typename property_t, typename next_receiver_t>
        struct  receiver
        {
            jst_t _jst;
            property_t _property;
            next_receiver_t _next_receiver;

        private:

            // we do not know where the observer is coming from?
            // and also when we move it through we might loose some information
            template <observable_stack subscriber_t>
            constexpr auto make_tree([[maybe_unused]] subscriber_t &subscriber)
            {
                if constexpr (libjst::is_resumable_v<std::remove_cvref_t<property_t>>)
                    return lazy_tree{libjst::root(_jst, libjst::window_size(_property)), subscriber};
                else
                    return lazy_tree{libjst::root(_jst, libjst::window_size(_property))};
            }

            template <traversable_journaled_sequence_tree haystack_t, typename algorithm_t>
            constexpr friend void tag_invoke(std::tag_t<libjst::set_value>,
                                             receiver && me,
                                             haystack_t && haystack,
                                             algorithm_t && algorithm) noexcept
            {
                // build the tree
                libjst::set_value(std::move(me._next_receiver),
                                  make_tree((haystack_t &&) haystack, libjst::window_size(algorithm)),
                                  (algorithm_t &&) algorithm);
            }

            template <typename cpo_t>
                requires std::invocable<cpo_t, next_receiver_t>
            constexpr friend void tag_invoke(cpo_t cpo, receiver && me) noexcept
            {
                cpo(std::move(me._next_receiver));
            }
        };

        template <typename predecessor_t, typename property_t>
        class searcher
        {
        private:
            predecessor_t _predecessor;
            property_t _property;

        public:

            // constructor
            template <typename _predecessor_t, typename _property_t>
            constexpr explicit searcher(_predecessor_t &&predecessor, _property_t &&property)
                noexcept(std::is_nothrow_constructible_v<predecessor_t, _predecessor_t> &&
                         std::is_nothrow_constructible_v<property_t, _property_t>) :
                _predecessor{(_predecessor_t &&)predecessor},
                _property{(_property_t &&)property}
            {}
        private:

            template <typename searcher_t, typename receiver_t>
            using operation_t = stateless_property_operation<jst::contrib::member_type_t<searcher_t, std::remove_cvref_t<predecessor_t>>,
                                                            jst::contrib::member_type_t<searcher_t, std::remove_cvref_t<property_t>>,
                                                            receiver_t>;

            template <typename searcher_t, typename receiver_t>
                requires (std::same_as<std::remove_cvref_t<searcher_t>, searcher>)
                        //   sender_to<std::remove_cvref_t<searcher_t>, receiver_t>)
            constexpr friend auto tag_invoke(std::tag_t<libjst::connect>, searcher_t && me, receiver_t && receiver)
                -> operation_t<searcher_t, receiver_t>
            {
                using fwd_jst_t = jst::contrib::member_type_t<searcher_t,     std::remove_cvref_t<predecessor_t>>;
                using fwd_property_t = jst::contrib::member_type_t<searcher_t, std::remove_cvref_t<property_t>>;
                return operation_t<searcher_t, receiver_t>{(fwd_jst_t &&)me._predecessor,
                                                            (fwd_property_t &&)me._property,
                                                            (receiver_t &&)receiver};
            }
        };

        inline constexpr struct _cpo
        {
            // we need to call it with some other parameter?
            template <typename searcher_t, typename property_t>
                // this must be a sender interface, how can we create this?
            constexpr auto operator()(searcher_t && searcher, property_t && property) const
                noexcept(std::is_nothrow_constructible_v<_searcher_on_lazy_tree::searcher<searcher_t, property_t>, searcher_t, property_t>)
                -> _searcher_on_lazy_tree::searcher<searcher_t, property_t>
            {
                return _searcher_on_lazy_tree::searcher<searcher_t, property_t>{(searcher_t &&)searcher,
                                                                                (property_t &&)property};
            }

            template <typename property_t>
            constexpr auto operator()(property_t && property) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, property_t>)
                -> jst::contrib::closure_result_t<_cpo, property_t>
            {
                return jst::contrib::make_closure(_cpo{}, (property_t &&)property);
            }
        } on_lazy_tree{};
    } // namespace _searcher_on_lazy_tree

    using _searcher_on_lazy_tree::on_lazy_tree;
} // namespace libjst
