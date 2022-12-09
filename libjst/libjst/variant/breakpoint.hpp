// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides breakpoint implementation.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <compare>
#include <concepts>
#include <string_view>

namespace libjst
{
    enum struct breakpoint_end : bool {
        right = 0,
        left = 1
    };

    class breakpoint {
    private:

        uint32_t _value      :31;
        uint32_t _end_marker :1;

    public:

        using value_type = uint32_t;

        constexpr breakpoint() noexcept : breakpoint{0u}
        {}

        constexpr breakpoint(uint32_t const value) noexcept : breakpoint{value, breakpoint_end::left}
        {}

        constexpr breakpoint(uint32_t const value, breakpoint_end const end_marker) noexcept :
            _value{value},
            _end_marker{static_cast<uint32_t>(end_marker)} // end marker == left => true => 1 > 0, so natural orering is right before left.
        {}

        constexpr uint32_t value() const noexcept {
            return _value;
        }

        constexpr breakpoint_end end_marker() const noexcept {
            return (is_left_end()) ? breakpoint_end::left : breakpoint_end::right;
        }

        constexpr bool is_left_end() const noexcept {
            return _end_marker == static_cast<uint32_t>(breakpoint_end::left);
        }

        constexpr bool is_right_end() const noexcept {
            return !is_left_end();
        }

        template <std::integral int_t>
        constexpr operator int_t() const noexcept {
            return static_cast<int_t>(value());
        }
    private:

        constexpr friend bool operator==(breakpoint const & lhs, breakpoint const & rhs) noexcept {
            return lhs <=> rhs == 0;
        }

        constexpr friend std::strong_ordering operator<=>(breakpoint const & lhs, breakpoint const & rhs) noexcept {
            if (auto cmp = lhs.value() <=> rhs.value(); cmp == 0) {
                return lhs.end_marker() <=> rhs.end_marker();
            } else {
                return cmp;
            }
        }
    };

    template <typename stream_t, typename breakpoint_t>
        requires std::same_as<std::remove_cvref_t<breakpoint_t>, breakpoint>
    inline stream_t & operator<<(stream_t & stream, breakpoint_t && bp) {
        using namespace std::literals;
        auto end_marker_string = [&] { return (bp.is_left_end()) ? "left end"sv : "right end"sv; };
        return stream << end_marker_string() << ": " << bp.value();
    }
}  // namespace libjst
