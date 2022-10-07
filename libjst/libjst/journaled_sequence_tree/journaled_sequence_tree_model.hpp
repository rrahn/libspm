// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides basic implementation of a journaled sequence tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <ranges>

#include <seqan3/range/concept.hpp>

#include <libjst/sequence_variant/concept.hpp>
#include <libjst/journaled_sequence_tree/serialiser_concept.hpp>

namespace libjst
{
    template <seqan3::sequence sequence_t, typename variant_store_t>
        requires std::ranges::random_access_range<sequence_t> &&
                 covered_sequence_variant<std::ranges::range_reference_t<variant_store_t>>
    class journaled_sequence_tree_model
    {
    private:

        using variant_value_t = std::ranges::range_value_t<variant_store_t>;

        std::reference_wrapper<sequence_t const> _base_sequence;
        variant_store_t _variant_store{};
        size_t _sequence_count{};
    public:

        explicit journaled_sequence_tree_model(sequence_t const & sequence, size_t const count) noexcept :
            _base_sequence{sequence},
            _sequence_count{count}
        {}

        bool insert(variant_value_t covered_variant)
        {  // initial contract
            assert(end_position(covered_variant) <= std::ranges::size(_base_sequence.get()));
            assert(std::ranges::size(libjst::coverage(covered_variant)) == _sequence_count);

            auto inserted = _variant_store.insert(std::move(covered_variant));
            return inserted != _variant_store.end();
        }

    private:

        template <typename archive_t>
        constexpr friend auto tag_invoke(std::tag_t<libjst::load>,
                                         journaled_sequence_tree_model & me,
                                         archive_t & archive)
        {
            libjst::load_extern(archive, me._base_sequence.get());
            archive(me._sequence_count, me._variant_store);
        }

        template <typename archive_t>
        constexpr friend auto tag_invoke(std::tag_t<libjst::save>,
                                         journaled_sequence_tree_model const & me,
                                         archive_t & archive)
        {
            libjst::save_extern(archive, me._base_sequence.get());
            archive(me._sequence_count, me._variant_store);
        }

        constexpr friend auto tag_invoke(std::tag_t<libjst::base_sequence>,
                                         journaled_sequence_tree_model const &me) noexcept -> sequence_t const &
        {
            return me._base_sequence.get();
        }

        constexpr friend auto tag_invoke(std::tag_t<libjst::variant_store>,
                                         journaled_sequence_tree_model const &me) noexcept -> variant_store_t const &
        {
            return me._variant_store;
        }

        constexpr friend size_t tag_invoke(std::tag_t<libjst::size>, journaled_sequence_tree_model const &me) noexcept
        {
            return me._sequence_count;
        }

        template <sequence_variant variant_t>
        constexpr auto end_position(variant_t && variant) const {
            return libjst::position(variant) + std::ranges::size(libjst::insertion(variant));
        }
    };
}  // namespace libjst
