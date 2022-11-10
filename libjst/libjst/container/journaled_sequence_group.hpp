// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::journaled_sequence_group
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>
#include <type_traits>

#include <seqan3/range/concept.hpp>

#include <libjst/variant/concept.hpp>

namespace libjst
{
    // what do we know about serialisation here if this is a view?
    template <std::ranges::view source_t, std::ranges::random_access_range variant_store_t>
        requires std::ranges::random_access_range<source_t> &&
                 seqan3::sequence<source_t> &&
                 std::ranges::sized_range<source_t> &&
                 covered_sequence_variant<std::ranges::range_reference_t<variant_store_t>>
    class journaled_sequence_group
    {
    private:

        using variant_value_t = std::ranges::range_value_t<variant_store_t>;

        source_t _source; // referenced source sequence
        variant_store_t _variant_store{}; // owned
        size_t _sequence_count{};
    public:

        journaled_sequence_group() = default; // depends on the source sequence wether it is default constructible.

        template <typename _source_t>
            requires (!std::same_as<std::remove_cvref_t<_source_t>, journaled_sequence_group> &&
                      std::ranges::viewable_range<_source_t>)
        explicit journaled_sequence_group(_source_t && source, size_t const count)
            noexcept(std::is_nothrow_constructible_v<std::views::all_t<_source_t>>) :
            _source{std::views::all((_source_t &&)source)},
            _sequence_count{count}
        {}

        template <typename _source_t>
            requires (!std::same_as<std::remove_cvref_t<_source_t>, journaled_sequence_group> &&
                      std::ranges::viewable_range<_source_t>)
        explicit journaled_sequence_group(_source_t &&source, variant_store_t variant_store)
            noexcept(std::is_nothrow_constructible_v<std::views::all_t<_source_t>> &&
                     std::is_nothrow_move_constructible_v<variant_store_t>) :
            _source{std::views::all((_source_t &&)source)},
            _variant_store{std::move(variant_store)}
        {
            _sequence_count = std::ranges::size(libjst::coverage(_variant_store[0]));
            for (auto && variant : _variant_store) {
                if (end_position(variant) > std::ranges::size(_source) ||
                    std::ranges::size(libjst::coverage(variant)) != size())
                    throw std::runtime_error{"Invalid variant store."};
            }
        }

        size_t size() const noexcept
        {
            return _sequence_count;
        }
    };

    template <std::ranges::viewable_range source_t, typename store_t>
    journaled_sequence_group(source_t &&, store_t) -> journaled_sequence_group<std::views::all_t<source_t>, store_t>;
}  // namespace libjst
