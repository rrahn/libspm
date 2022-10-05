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

#include <ranges>
#include <vector>

#include <libjst/journaled_sequence_tree/concept.hpp>
#include <libjst/sequence_variant/concept.hpp>

namespace libio
{
    // might define my variant wrapper here!
    template <libjst::covered_sequence_variant variant_t>
    class reverse_variant
    {
        using position_t = variant_position_t<variant_t>;

        variant_t const & _wrappee;
        position_t _position;

    public:
        reverse_variant(variant_t const & variant, size_t const base_sequence_size) :
            _wrappee{variant}, _position{base_sequence_size - libjst::position(variant)}
        {}

    private:

        template <typename cpo_t>
            requires std::invocable<cpo_t, variant_t const &>;
        constexpr friend auto tag_invoke(cpo_t cpo, reverse_variant const &me)
            noexcept(std::is_nothrow_invocable_v<cpo_t, variant_t const &>)
            -> std::invoke_result_t<cpo_t, variant_t const &>
        {
            return cpo(me._wrappee);
        }

        constexpr friend auto tag_invoke(std::tag_t<libjst::position>, reverse_variant const &me)
            noexcept(std::is_nothrow_invocable_v<std::tag_t<libjst::position>, variant_t const &>)
            -> std::invoke_result_t<std::tag_t<libjst::position>, variant_t const &>
        {
            return me._position;
        }

        constexpr friend auto tag_invoke(std::tag_t<libjst::insertion>, reverse_variant const &me)
            noexcept(std::is_nothrow_invocable_v<std::tag_t<libjst::insertion>, variant_t const &>)
            -> decltype(std::declval<std::invoke_result_t<libjst::insertion, variant_t const &>>() | std::views::reverse)
        {
            return libjst::insertion(me._wrappee) | std::views::reverse;
        }
    };

    template <typename jst_t>
    class journaled_sequence_tree_backward : traversable_jst_base
    {
    private:
        using variant_store_t = variant_store_t<jst_t>;
        using variant_reference_t = std::ranges::range_reference_t<variant_store_t>;
        using reverse_variant_t = reverse_variant<variant_reference_t>;
        using sorted_variant_store_t = std::vector<reverse_variant_t>;

        [[no_unique_address]] jst_t const &_jst;
        sorted_variant_store_t _event_queue;

    public:

        explicit journaled_sequence_tree_backward(jst_t const & jst) : _jst{jst}
        {
            _event_queue.resize(libjst::variant_store(_jst).size());
            std::ranges::copy(libjst::variant_store(_jst) | std::views::reverse
                                                          | std::views::transform([&] (variant_reference_t variant)
            {
                return reverse_variant_t{variant, std::ranges::size(libjst::base_sequence(_jst))};
            }), std::ranges::begin(_event_queue));

            // this is now the same for the ordering.
            auto effective_size = [] (auto const &variant) {
                return std::ranges::size(lib::insertion(variant)) - lib::deletion(variant);
            };

            auto cmp = [&] (auto const &lhs, auto const &rhs)
            {
                // we need another variant adator proxy!
                return (libjst::position(lhs) < libjst::position(rhs)) ||
                       (libjst::position(lhs) == libjst::position(rhs) && effective_size(lhs) > effective_size(rhs));
            };

            std::ranges::stable_sort(_event_queue, cmp);
            // assert(std::ranges::is_sorted(_event_queue)); // make sure the range is sorted.
        }

    private:
        template <typename cpo_t>
            requires std::invocable<cpo_t, jst_t const &>;
        constexpr friend auto tag_invoke(cpo_t cpo, journaled_sequence_tree_backward const &me)
            noexcept(std::is_nothrow_invocable_v<cpo_t, jst_t const &>)
            -> std::invoke_result_t<cpo_t, jst_t const &>
        {
            return cpo(me._jst);
        }

        constexpr friend auto tag_invoke(std::tag_t<libjst::base_sequence>,
                                         journaled_sequence_tree_backward const &me) noexcept -> variant_store_t const &
        {
            return std::views::reverse(libjst::base_sequence(me._jst));
        }

        constexpr friend auto tag_invoke(std::tag_t<libjst::variant_store>,
                                         journaled_sequence_tree_backward const &me) noexcept -> variant_store_t const &
        {
            return me._event_queue;
        }
    };

}  // namespace libio
