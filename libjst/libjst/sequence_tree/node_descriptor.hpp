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

    enum class node_descriptor_id {
        nil = 0,
        reference = 1,
        alternate = 2,
        on_alternate_path = 4,
        // branching = 0,
        // reference node categories
        first_left = 8,
        first_right = 16,
        second_left = 32,
        second_right = 64,
        second_first_right = 128,
        first_breakpoint_mask = first_left | first_right,
        second_breakpoint_mask = second_left | second_right | second_first_right
    };

    enum class node_state {
        A,
        B,
        C,
        D,
        E,
        F,
        G,
        H,
        I,
        final
    };

} // namespace libjst

namespace seqan3 {

    template <>
    constexpr bool add_enum_bitwise_operators<libjst::node_descriptor_id> = true;
} // namespace seqan3

namespace libjst {

    class node_descriptor {
    private:
        node_descriptor_id _value{};

        static constexpr node_descriptor_id unset_branching_mask{[] {
            using seqan3::operator|;
            return node_descriptor_id::reference | node_descriptor_id::alternate | node_descriptor_id::on_alternate_path;
        } ()};

    public:

        constexpr node_descriptor() = default;

        constexpr bool from_reference() const noexcept {
            using seqan3::operator&;
            return (_value & node_descriptor_id::reference) == node_descriptor_id::reference;
        }

        constexpr bool from_alternate() const noexcept {
            using seqan3::operator&;
            return (_value & node_descriptor_id::alternate) == node_descriptor_id::alternate;
        }

        constexpr bool is_branching() const noexcept {
            return get_second_breakpoint_id() == node_descriptor_id::second_left;
        }

        constexpr bool on_alternate_path() const noexcept {
            using seqan3::operator&;
            return (_value & node_descriptor_id::on_alternate_path) == node_descriptor_id::on_alternate_path;
        }

        constexpr void set_reference() noexcept {
            using seqan3::operator|;
            using seqan3::operator&;
            _value = node_descriptor_id::reference | (_value & node_descriptor_id::on_alternate_path);
        }

        constexpr void set_alternate() noexcept {
            using seqan3::operator|;
            _value = node_descriptor_id::alternate | node_descriptor_id::on_alternate_path;
        }

        // constexpr void set_branching() noexcept {
        //     using seqan3::operator|=;
        //     _value |= node_descriptor_id::branching;
        // }

        // constexpr void unset_branching() noexcept {
        //     using seqan3::operator&=;
        //     _value &= unset_branching_mask;
        // }

        constexpr node_descriptor_id get_first_breakpoint_id() const noexcept {
            using seqan3::operator&;
            return _value & node_descriptor_id::first_breakpoint_mask;
        }

        constexpr node_descriptor_id get_second_breakpoint_id() const noexcept {
            using seqan3::operator&;
            return _value & node_descriptor_id::second_breakpoint_mask;
        }

        constexpr void set_first_breakpoint_id(node_descriptor_id const bp_id) noexcept {
            using seqan3::operator&=;
            using seqan3::operator|=;
            using seqan3::operator&;
            using seqan3::operator~;

            assert((bp_id & node_descriptor_id::first_breakpoint_mask) != node_descriptor_id::nil);
            _value &= ~node_descriptor_id::first_breakpoint_mask;
            _value |= bp_id;
            // _value |= (bp_id & node_descriptor_id::first_breakpoint_mask);
        }

        constexpr void set_second_breakpoint_id(node_descriptor_id const bp_id) noexcept {
            using seqan3::operator&=;
            using seqan3::operator|=;
            using seqan3::operator&;
            using seqan3::operator~;

            assert((bp_id & node_descriptor_id::second_breakpoint_mask) != node_descriptor_id::nil);
            _value &= ~node_descriptor_id::second_breakpoint_mask;
            _value |= bp_id;
        }

    private:
        constexpr friend bool operator==(node_descriptor const & lhs, node_descriptor const & rhs) noexcept = default;
    };
}  // namespace libjst
