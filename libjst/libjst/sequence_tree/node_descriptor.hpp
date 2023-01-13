// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides node descriptor.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/add_enum_bitwise_operators.hpp>

namespace libjst
{
    enum class breakpoint_state : uint8_t {
        nil              = 0b0000,
        left_begin       = 0b1000,
        left_end         = 0b0100,
        right_begin      = 0b0010,
        right_end        = 0b0001,
    };

    enum class node_state : uint32_t {
        nil = 0b00000,
        last = 0b10000,
/*A*/   branching_after_left_end        = static_cast<uint32_t>(breakpoint_state::left_end) |
                                          static_cast<uint32_t>(breakpoint_state::right_begin),
/*B*/   last_branching_after_left_end   = static_cast<uint32_t>(breakpoint_state::left_end) |
                                          static_cast<uint32_t>(breakpoint_state::right_begin) | last,
/*E*/   branching_after_left_begin      = static_cast<uint32_t>(breakpoint_state::left_begin) |
                                          static_cast<uint32_t>(breakpoint_state::right_begin),
/*F*/   last_branching_after_left_begin = static_cast<uint32_t>(breakpoint_state::left_begin) |
                                          static_cast<uint32_t>(breakpoint_state::right_begin) | last,
/*D*/   non_branching_left_only         = static_cast<uint32_t>(breakpoint_state::left_begin) |
                                          static_cast<uint32_t>(breakpoint_state::left_end),
/*C*/   last_non_branching_left_only    = static_cast<uint32_t>(breakpoint_state::left_begin) |
                                          static_cast<uint32_t>(breakpoint_state::left_end) | last,
/*H*/   non_branching_including_left    = static_cast<uint32_t>(breakpoint_state::left_begin) |
                                          static_cast<uint32_t>(breakpoint_state::right_end),
/*G*/   non_branching_after_left        = static_cast<uint32_t>(breakpoint_state::left_end) |
                                          static_cast<uint32_t>(breakpoint_state::right_end),
        variant                         = static_cast<uint32_t>(breakpoint_state::left_begin) |
                                          static_cast<uint32_t>(breakpoint_state::left_end) |
                                          static_cast<uint32_t>(breakpoint_state::right_end),
    };

} // namespace libjst

namespace seqan3 {
    template <>
    constexpr bool add_enum_bitwise_operators<libjst::breakpoint_state> = true;

    template <>
    constexpr bool add_enum_bitwise_operators<libjst::node_state> = true;
} // namespace seqan3

namespace libjst {

    class breakpoint_descriptor {
    private:
        breakpoint_state _state{};
    public:

        constexpr breakpoint_descriptor() = default;
        constexpr explicit breakpoint_descriptor(breakpoint_state state) noexcept : _state{state}
        {}

        constexpr breakpoint_descriptor& operator=(breakpoint_state state) noexcept {
            _state = state;
            return *this;
        }

        constexpr explicit operator breakpoint_state() const noexcept {
            return _state;
        }

        constexpr bool from_left_begin() const noexcept {
            using seqan3::operator&;
            return (_state & breakpoint_state::left_begin) == breakpoint_state::left_begin;
        }

        constexpr bool from_left_end() const noexcept {
            using seqan3::operator&;
            return (_state & breakpoint_state::left_end) == breakpoint_state::left_end;
        }

        constexpr bool from_right_begin() const noexcept {
            using seqan3::operator&;
            return (_state & breakpoint_state::right_begin) == breakpoint_state::right_begin;
        }

        constexpr bool from_right_end() const noexcept {
            using seqan3::operator&;
            return (_state & breakpoint_state::right_end) == breakpoint_state::right_end;
        }
    private:

        constexpr friend bool operator==(breakpoint_descriptor const &, breakpoint_descriptor const &) noexcept = default;
    };

    class node_descriptor {
    private:
        breakpoint_descriptor _left_break{};
        breakpoint_descriptor _right_break{};
        bool _last{};
        bool _on_alternate_path{};

    public:

        constexpr node_descriptor() = default;
        constexpr explicit node_descriptor(node_state const from_state) noexcept {
            activate_state(from_state);
        }

        constexpr node_descriptor & operator=(node_state const from_state) noexcept {
            activate_state(from_state);
            return *this;
        }

        constexpr void activate_state(node_state const from_state) noexcept {
            set_left_break(from_state);
            set_right_break(from_state);
            set_last_state(from_state);
        }

        constexpr bool from_reference() const noexcept {
            return !from_variant();
        }

        constexpr bool from_variant() const noexcept {
            return _left_break.from_left_begin() && _left_break.from_left_end() && _right_break.from_right_end();
        }

        constexpr bool is_branching() const noexcept {
            return _right_break.from_right_begin();
        }

        constexpr bool on_alternate_path() const noexcept {
            return _on_alternate_path;
        }

        // not sure what to do here! actually not needed.
        constexpr void set_reference(breakpoint_state left_state,
                                     breakpoint_state right_state,
                                     bool last_state = false) noexcept {
            _left_break = left_state;
            _right_break = right_state;
            _last = last_state;
            assert(from_reference());
        }

        constexpr void set_variant() noexcept {
            using seqan3::operator|;
            _left_break = breakpoint_state::left_begin | breakpoint_state::left_end | breakpoint_state::right_end;
            _right_break = _left_break;
            _on_alternate_path = true;
        }

        constexpr breakpoint_descriptor const & left_break() const noexcept {
            return _left_break;
        }

        constexpr breakpoint_descriptor const & right_break() const noexcept {
            return _right_break;
        }

        constexpr operator node_state() const noexcept {
            using seqan3::operator|;
            return node_state{static_cast<uint32_t>(static_cast<breakpoint_state>(left_break()) |
                                                    static_cast<breakpoint_state>(right_break()))} |
                   ((_last) ? node_state::last : node_state::nil);
        }

    private:

        constexpr void set_left_break(node_state const from_state) noexcept {
            using seqan3::operator|;
            switch (from_state) {
                case node_state::branching_after_left_end: [[fallthrough]];
                case node_state::last_branching_after_left_end: [[fallthrough]];
                case node_state::non_branching_after_left: _left_break = breakpoint_state::left_end; break;
                case node_state::branching_after_left_begin: [[fallthrough]];
                case node_state::last_branching_after_left_begin: [[fallthrough]];
                case node_state::non_branching_left_only: [[fallthrough]];
                case node_state::last_non_branching_left_only: [[fallthrough]];
                case node_state::non_branching_including_left: _left_break = breakpoint_state::left_begin; break;
                case node_state::variant: set_variant(); break;
                default: _left_break = breakpoint_state::nil;
            }
        }

        constexpr void set_right_break(node_state const from_state) noexcept {
            using seqan3::operator|;
            switch (from_state) {
                case node_state::branching_after_left_end: [[fallthrough]];
                case node_state::last_branching_after_left_end: [[fallthrough]];
                case node_state::branching_after_left_begin: [[fallthrough]];
                case node_state::last_branching_after_left_begin: _right_break = breakpoint_state::right_begin; break;
                case node_state::non_branching_left_only: [[fallthrough]];
                case node_state::last_non_branching_left_only: _right_break = breakpoint_state::left_end; break;
                case node_state::non_branching_after_left: [[fallthrough]];
                case node_state::non_branching_including_left: _right_break = breakpoint_state::right_end; break;
                case node_state::variant: set_variant(); break;
                default: _right_break = breakpoint_state::nil;
            }
        }

        constexpr void set_last_state(node_state const from_state) noexcept {
            switch (from_state) {
                case node_state::last_branching_after_left_end: [[fallthrough]];
                case node_state::last_branching_after_left_begin: [[fallthrough]];
                case node_state::last_non_branching_left_only: _last = true; break;
                default: _last = false;
            }
        }

        constexpr friend bool operator==(node_descriptor const & lhs, node_descriptor const & rhs) noexcept = default;
    };
}  // namespace libjst
