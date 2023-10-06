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

#include <libcontrib/type_traits.hpp>

#include <libjst/variant/concept.hpp>

namespace libjst
{
    class breakpoint {
    private:

        uint32_t _value      :31;
        uint32_t _end_marker :1;

        uint32_t _low{};
        uint32_t _high{};

    public:

        using value_type = uint32_t;

        constexpr breakpoint() noexcept : breakpoint{0u, breakpoint_end::left}
        {}
        constexpr breakpoint(breakpoint const & other) noexcept :
            _value{other._value},
            _end_marker{other._end_marker},
            _low{other._low},
            _high{other._high}
        {}
        constexpr breakpoint(breakpoint && other) noexcept : breakpoint{std::as_const(other)}
        {}
        constexpr breakpoint& operator=(breakpoint other) noexcept
        {
            _value = other._value;
            _end_marker = other._end_marker;
            _low = other._low;
            _high = other._high;
            return *this;
        }
        ~breakpoint() = default;

        constexpr breakpoint(uint32_t const value) noexcept : breakpoint{value, breakpoint_end::left}
        {}

        constexpr breakpoint(uint32_t const value, breakpoint_end const end_marker) noexcept :
            _value{value},
            _end_marker{static_cast<uint32_t>(end_marker)} // end marker == left => true => 1 > 0, so natural orering is right before left.
        {}

        constexpr breakpoint(value_type const low, std::size_t const count) noexcept :
            _value{low},
            _end_marker{1},
            _low{low},
            _high{low + static_cast<value_type>(count)}
        {
            assert(static_cast<std::size_t>(low) + count <= std::numeric_limits<value_type>::max());
        }

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
        explicit constexpr operator int_t() const noexcept {
            return static_cast<int_t>(value());
        }

        template <typename output_archive_t>
        uint32_t save_minimal(output_archive_t const &) const
        {
            uint32_t tmp = _end_marker;
            tmp = (tmp << 31) | _value; // make space for the position.
            return tmp;
        }

        template <typename input_archive_t>
        void load_minimal(input_archive_t const &, uint32_t const &tmp)
        {
            _value = tmp;
            _end_marker = (tmp >> 31);
        }

    private:

        template <typename me_t>
            requires std::same_as<std::remove_cvref_t<me_t>, breakpoint>
        constexpr friend auto tag_invoke(libjst::tag_t<libjst::low_breakend>, me_t && me) noexcept
            -> jst::contrib::member_type_t<me_t, value_type>
        {
            return static_cast<jst::contrib::member_type_t<me_t, value_type>>(me._low);
        }

        template <typename me_t>
            requires std::same_as<std::remove_cvref_t<me_t>, breakpoint>
        constexpr friend auto tag_invoke(libjst::tag_t<libjst::high_breakend>, me_t && me) noexcept
            -> jst::contrib::member_type_t<me_t, value_type>
        {
            return static_cast<jst::contrib::member_type_t<me_t, value_type>>(me._high);
        }

        template <typename me_t>
            requires std::same_as<std::remove_cvref_t<me_t>, breakpoint>
        constexpr friend auto tag_invoke(libjst::tag_t<libjst::breakpoint_span>, me_t && me) noexcept
            -> value_type
        {
            return libjst::high_breakend(me) - libjst::low_breakend(me);
        }

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
        return stream << "[" << libjst::low_breakend(bp) << ".." << libjst::high_breakend(bp) << ")";
    }
}  // namespace libjst
