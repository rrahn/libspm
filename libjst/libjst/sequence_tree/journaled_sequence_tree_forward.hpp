// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of the journaled sequence tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/closure_adaptor.hpp>
#include <libcontrib/copyable_box.hpp>

#include <seqan3/utility/type_traits/lazy_conditional.hpp>
#include <seqan3/range/views/type_reduce.hpp>

#include <libjst/journaled_sequence_tree/concept.hpp>
#include <libjst/journaled_sequence_tree/serialiser_concept.hpp>
#include <libjst/journaled_sequence_tree/serialiser_delegate.hpp>
#include <libjst/sequence_variant/variant_store_sorted.hpp>
#include <libjst/traversal/jst_node_base.hpp>
#include <libjst/traversal/jst_node_value.hpp>

namespace libjst
{

    namespace detail
    {
        template <typename jst_t>
        using maybe_reference_wrapper_t = std::conditional_t<std::is_lvalue_reference_v<jst_t>,
                                                             std::reference_wrapper<std::remove_reference_t<jst_t>>,
                                                             jst_t>;
        template <typename jst_t>
        class jst_box : public jst::contrib::copyable_box<maybe_reference_wrapper_t<jst_t>>
        {
            using base_t = jst::contrib::copyable_box<maybe_reference_wrapper_t<jst_t>>;
            // TODO: add serialiser functionality
        public:

            using base_t::base_t;

            constexpr jst_t & operator*() noexcept(noexcept(*std::declval<base_t &>()))
            {
                return *static_cast<base_t &>(*this);
            }

            constexpr jst_t const & operator*() const noexcept(noexcept(*std::declval<base_t const &>()))
            {
                return *static_cast<base_t const &>(*this);
            }
        };

    } // namespace detail

    template <typename jst_t>
    class journaled_sequence_tree_forward_ : private traversable_jst_base
    {
    private:
        using store_t = variant_store_t<jst_t>;
        using sorted_store_t = variant_store_sorted<store_t const>;

        class ordered_variant;

        // we need to put it into a moveable-box!
        // then we do not need to store a reference object here.
        [[no_unique_address]] detail::jst_box<jst_t> _jst{};
        sorted_store_t _event_queue{};

    public:

        journaled_sequence_tree_forward_() = default;
        journaled_sequence_tree_forward_(journaled_sequence_tree_forward_ const &) = default;
        journaled_sequence_tree_forward_(journaled_sequence_tree_forward_ &&) = default;
        journaled_sequence_tree_forward_ & operator=(journaled_sequence_tree_forward_ const &) = default;
        journaled_sequence_tree_forward_ & operator=(journaled_sequence_tree_forward_ &&) = default;

        template <typename _jst_t>
            requires (!std::same_as<std::remove_cvref_t<_jst_t>, journaled_sequence_tree_forward_> &&
                       std::same_as<std::remove_cvref_t<_jst_t>, std::remove_cvref_t<jst_t>>)
        explicit journaled_sequence_tree_forward_(_jst_t &&jst) :
            _jst{(_jst_t &&) jst},
            _event_queue{libjst::variant_store(*_jst)}
        {}

    private:
        template <typename archive_t>
        constexpr friend auto tag_invoke(std::tag_t<libjst::load>,
                                         journaled_sequence_tree_forward_ & me,
                                         archive_t & archive)
        {
            libjst::load_extern(archive, *me._jst);
            archive(me._event_queue);
        }

        template <typename archive_t>
        constexpr friend auto tag_invoke(std::tag_t<libjst::save>,
                                         journaled_sequence_tree_forward_ const & me,
                                         archive_t & archive)
        {
            libjst::save_extern(archive, *me._jst);
            archive(me._event_queue);
        }

        template <typename cpo_t>
            requires std::invocable<cpo_t, jst_t const &>
        constexpr friend auto tag_invoke(cpo_t cpo, journaled_sequence_tree_forward_ const &me)
            noexcept(std::is_nothrow_invocable_v<cpo_t, jst_t const &>)
            -> std::invoke_result_t<cpo_t, jst_t const &>
        {
            return cpo(*me._jst);
        }

        constexpr friend auto tag_invoke(std::tag_t<libjst::variant_store>,
                                         journaled_sequence_tree_forward_ const &me) noexcept
            -> sorted_store_t const &
        {
            return me._event_queue;
        }

        // Very complex type definitions.
        using base_reduce_t = decltype(std::declval<base_sequence_t<jst_t const &> const &>() | seqan3::views::type_reduce);
        using variant_t = std::ranges::range_value_t<sorted_store_t>;
        using position_t = variant_position_t<variant_t>;
        using coverage_t = variant_coverage_t<variant_t>;
        using ins_reduce_t = decltype(std::declval<variant_insertion_t<variant_t const &>>() | seqan3::views::type_reduce);

        using entry_sequence_t = std::common_type_t<base_reduce_t, ins_reduce_t>;
        using journal_t = journal<position_t, entry_sequence_t>;
        using node_value_t = jst_node_value<journal_t, coverage_t>;
        using node_t = jst_node_base<node_value_t, std::ranges::iterator_t<sorted_store_t const>>;

        template <typename branch_state_t>
        constexpr friend auto tag_invoke(std::tag_t<libjst::root>,
                                         journaled_sequence_tree_forward_ const &me,
                                         branch_state_t && initial_branch_state,
                                         size_t const window_size) noexcept -> node_t
        {
            // we might need to wrap this stuff anyway!
            auto base_view = libjst::base_sequence(*me._jst) | seqan3::views::type_reduce;

            node_value_t value{base_view, coverage_t(true, libjst::size(*me._jst))};
            return node_t{std::move(value),
                          std::ranges::begin(me._event_queue),
                          std::ranges::end(me._event_queue),
                          window_size};
        }

    };

    template <typename jst_t>
    journaled_sequence_tree_forward_(jst_t &&) -> journaled_sequence_tree_forward_<jst_t>;

    namespace _forward_jst
    {
        inline constexpr struct _fn
        {
            template <journaled_sequence_tree_c jst_model_t>
            constexpr auto operator()(jst_model_t &&jst) const
                -> journaled_sequence_tree_forward_<jst_model_t>
            {
                return journaled_sequence_tree_forward_<jst_model_t>{(jst_model_t &&)jst};
            }

            template <typename ...args_t>
                requires (sizeof...(args_t) == 0)
            constexpr auto operator()(args_t &&...) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, _fn, args_t...>)
                -> jst::contrib::closure_result_t<_fn, args_t...>
            {
                return jst::contrib::make_closure(_fn{});
            }

        } forward_jst;
    } // namespace _forward_jst

    using _forward_jst::forward_jst;
}  // namespace libjst
