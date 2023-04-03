// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides common sequence interface for different sequence variant encodings.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>
#include <span>
#include <variant>

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/utility/detail/multi_invocable.hpp>

namespace libjst
{
    template <std::ranges::contiguous_range source_t>
    class delta_sequence_variant : public std::ranges::view_base {

        using value_type = std::ranges::range_value_t<source_t>;
        using deletion_sequence_type = std::array<value_type, 0>;
        using snv_sequence_type = std::array<value_type, 1>;
        using span_type = std::span<value_type const>;

        using sequence_type = std::variant<deletion_sequence_type, snv_sequence_type, span_type>;

        sequence_type _sequence{};
    public:

        using iterator = std::ranges::iterator_t<span_type>;

        constexpr delta_sequence_variant() noexcept = default;

        explicit constexpr delta_sequence_variant(value_type snv) noexcept : _sequence{snv_sequence_type{snv}}
        {}

        explicit constexpr delta_sequence_variant(source_t const & insertion) noexcept : _sequence{span_type{insertion}}
        {}

        constexpr iterator begin() const noexcept {
            return std::ranges::begin(get_active_span());
        }

        constexpr iterator end() const noexcept {
            return std::ranges::end(get_active_span());
        }

    private:

        constexpr span_type get_active_span() const noexcept {
            return std::visit(seqan3::detail::multi_invocable{
                [] (span_type insertion_sequence) { return insertion_sequence; },
                [] (auto const & array_sequence) { return span_type{array_sequence}; }
            }, _sequence);
        }
    };

}  // namespace libjst
