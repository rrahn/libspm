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

#include <algorithm>
#include <ranges>
#include <vector>

#include <cereal/types/vector.hpp>

#include <seqan3/range/views/type_reduce.hpp>

#include <libjst/sequence_variant/variant_store_iterator.hpp>
#include <libjst/journaled_sequence_tree/concept.hpp>
#include <libjst/journaled_sequence_tree/serialiser_concept.hpp>

namespace libjst
{

    template <typename jst_t>
    class journaled_sequence_tree_forward_ : private traversable_jst_base
    {
    private:
        using store_t = variant_store_t<jst_t const &>;
        using variant_reference_t = std::ranges::range_reference_t<store_t>;
        using position_t = variant_position_t<variant_reference_t>;

        class ordered_variant;

        [[no_unique_address]] jst_t const *_jst{};
        std::vector<position_t> _event_queue{};

    public:

        journaled_sequence_tree_forward_() = default;
        // journaled_sequence_tree_forward_(journaled_sequence_tree_forward_ const &) = default;
        // journaled_sequence_tree_forward_(journaled_sequence_tree_forward_ &&) = default;
        // journaled_sequence_tree_forward_ &operator=(journaled_sequence_tree_forward_ const &) = default;
        // journaled_sequence_tree_forward_ &operator=(journaled_sequence_tree_forward_ &&) = default;
        explicit journaled_sequence_tree_forward_(jst_t const & jst) : _jst{std::addressof(jst)}
        {
            position_t store_size = std::ranges::size(libjst::variant_store(jst));
            _event_queue.resize(store_size);
            std::ranges::copy(std::views::iota(0u, store_size), std::ranges::begin(_event_queue));

            auto effective_size = [] (auto const &variant) {
                return std::ranges::size(libjst::insertion(variant)) - libjst::deletion(variant);
            };

            auto cmp = [&] (position_t const &left_position, position_t const &right_position) -> bool
            {
                variant_reference_t lhs = libjst::variant_store(jst)[left_position];
                variant_reference_t rhs = libjst::variant_store(jst)[right_position];
                return (libjst::position(lhs) < libjst::position(rhs)) ||
                       (libjst::position(lhs) == libjst::position(rhs) && effective_size(lhs) > effective_size(rhs));
            };

            std::ranges::stable_sort(_event_queue.begin(), _event_queue.end(), cmp); // now sort based on the indices.
            // assert(std::ranges::is_sorted(_event_queue)); // make sure the range is sorted.
        }

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
            assert(me._jst != nullptr);
            return cpo(*me._jst);
        }

        constexpr friend decltype(auto) tag_invoke(std::tag_t<libjst::variant_store>,
                                                   journaled_sequence_tree_forward_ const &me) noexcept
        {
            assert(me._jst != nullptr);
            return me._event_queue | std::views::transform([&] (position_t const position) {
                return ordered_variant{libjst::variant_store(*me._jst), position};
            });
        }

        constexpr friend auto tag_invoke(std::tag_t<libjst::base_sequence>,
                                         journaled_sequence_tree_forward_ const &me) noexcept
            -> decltype(std::declval<base_sequence_t<jst_t> const &>() | seqan3::views::type_reduce)
        {
            assert(me._jst != nullptr);
            return libjst::base_sequence(*me._jst) | seqan3::views::type_reduce;
        }
    };

    template <typename jst_t>
    class journaled_sequence_tree_forward_<jst_t>::ordered_variant
    {
    private:
        store_t const & _store{};
        position_t _position{};

    public:

        ordered_variant() = delete;
        explicit ordered_variant(store_t const & store, position_t const position) noexcept :
            _store{store},
            _position{position}
        {
        }

    private:

        template <typename cpo_t>
            requires (std::tag_invocable<cpo_t, std::ranges::range_reference_t<store_t const>>)
        constexpr friend auto tag_invoke(cpo_t cpo, ordered_variant const &me)
            -> std::tag_invoke_result_t<cpo_t, std::ranges::range_reference_t<store_t const>>
        {
            return std::tag_invoke(cpo, me._store[me._position]);
        }
    };

    // semantic:
    //  - no ambiguous covered variants
    //  - this is now a traversable jst => meaning that we require that its elements are ordered.
    //  - ser
    //

    // I can now independently get the objects from the values
    // now I could say search -> that stays as long inside the operation until done.
    // journaled_seuence_tree_forward.search(pattern_t, delegate?) // what can we do now?
    // the search and the operation is not a state!
    // we need to extract the search pattern
    // so another step is the handling of the search -> what can we do to make it searchable?
    // we have another algorithm
    // related classes such as the journaled sequence trees
    // interesting
}  // namespace libjst
