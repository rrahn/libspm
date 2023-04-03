// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a packed key type for the breakend dictionary of an rcms object.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <type_traits>

#include <seqan3/core/concept/cereal.hpp>

namespace libjst
{
    enum class indel_breakend_kind {
        deletion_high = 0b100,
        insertion_low = 0b101,
        deletion_low = 0b110
    };

    template <typename position_t = uint32_t>
    class packed_breakend_key {
    private:

        static constexpr position_t indel_mask{0b100};
        static constexpr position_t snv_mask{0b011};

        position_t _code : 3;
        position_t _position : 29;
    public:

        using underlying_type = position_t;

        constexpr packed_breakend_key() noexcept : _code{0}, _position{0}
        {}

        constexpr packed_breakend_key(indel_breakend_kind indel_kind, underlying_type position) noexcept :
            _code{static_cast<underlying_type>(indel_kind)},
            _position{position}
        {}

        constexpr packed_breakend_key(uint8_t snv_value, underlying_type position) noexcept :
            _code{static_cast<underlying_type>(snv_value) & snv_mask},
            _position{position}
        {}

        constexpr bool is_indel() const noexcept {
            return _code & indel_mask;
        }

        constexpr indel_breakend_kind indel_kind() const noexcept {
            assert(is_indel());
            return static_cast<indel_breakend_kind>(_code);
        }

        constexpr underlying_type snv_value() const noexcept {
            assert(!is_indel());
            return _code;
        }

        constexpr underlying_type position() const noexcept {
            return _position;
        }

        template <typename visitor_t>
        constexpr auto visit(visitor_t && fn) const
            noexcept(std::is_nothrow_invocable_v<visitor_t &&, indel_breakend_kind> &&
                     std::is_nothrow_invocable_v<visitor_t &&, underlying_type>)
            -> std::common_type_t<std::invoke_result_t<visitor_t &&, indel_breakend_kind>,
                                  std::invoke_result_t<visitor_t &&, underlying_type>>
        {
            if (is_indel()) {
                return std::invoke((visitor_t &&)fn, indel_kind());
            } else {
                return std::invoke((visitor_t &&)fn, snv_value());
            }
        }

        // ----------------------------------------------------------------------------
        // Serialisation
        // ----------------------------------------------------------------------------

        template <seqan3::cereal_input_archive archive_t>
        void load(archive_t & iarchive)
        {
            position_t packed_key{};
            iarchive(packed_key);
            _code = packed_key >> ((sizeof(position_t) * 8) - 3ull);
            _position = packed_key;
        }

        template <seqan3::cereal_output_archive archive_t>
        void save(archive_t & oarchive) const
        {
            position_t packed_key = _code << ((sizeof(position_t) * 8) - 3ull);
            packed_key |= _position;
            oarchive(packed_key);
        }


    private:

        constexpr underlying_type code() const noexcept {
            return _code;
        }

        constexpr friend bool operator==(packed_breakend_key const &, packed_breakend_key const &) = default;

        constexpr friend std::strong_ordering operator<=>(packed_breakend_key const & lhs,
                                                          packed_breakend_key const & rhs) noexcept {
            if (auto position_cmp = lhs.position() <=> rhs.position(); position_cmp == 0) {
                if (lhs.is_indel() == rhs.is_indel()) {
                    return lhs.code() <=> rhs.code();
                } else {
                    if (lhs.is_indel()) {
                        assert(!rhs.is_indel());
                        return (lhs.indel_kind() == indel_breakend_kind::deletion_low) ?
                                    std::strong_ordering::greater :
                                    std::strong_ordering::less;
                    } else {
                        assert(!lhs.is_indel());
                        return (rhs.indel_kind() == indel_breakend_kind::deletion_low) ?
                                    std::strong_ordering::less :
                                    std::strong_ordering::greater;
                    }
                }
            } else {
                return position_cmp;
            }
        }
    };

}  // namespace libjst

namespace std {
    template <typename position_t>
    struct hash<libjst::packed_breakend_key<position_t>> {
        constexpr std::size_t operator()(libjst::packed_breakend_key<position_t> const & key) const noexcept {
            return static_cast<std::size_t>(key.position()) |
                    key.visit([&] (auto const & code) { return static_cast<std::size_t>(code) << 29; });
        }
    };
}
