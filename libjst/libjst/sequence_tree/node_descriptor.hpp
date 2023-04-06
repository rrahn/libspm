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

#include <iosfwd>
#include <string_view>

#include <seqan3/core/add_enum_bitwise_operators.hpp>
#include <seqan3/core/concept/cereal.hpp>

namespace libjst
{
    enum class node_state : uint8_t {
        nil         = 0b00000,
        left_begin  = 0b10000,
        left_end    = 0b01000,
        right_begin = 0b00100,
        right_end   = 0b00010,
        last        = 0b00001,
/*A*/   branching_after_left_end        = left_end   | right_begin,
/*B*/   last_branching_after_left_end   = left_end   | right_begin | last,
/*E*/   branching_after_left_begin      = left_begin | right_begin,
/*F*/   last_branching_after_left_begin = left_begin | right_begin | last,
/*D*/   non_branching_left_only         = left_begin | left_end,
/*C*/   last_non_branching_left_only    = left_begin | left_end    | last,
/*H*/   non_branching_including_left    = left_begin | right_end,
/*G*/   non_branching_after_left        = left_end   | right_end,
        variant                         = left_begin | left_end    | right_end,
        reference                       = nil,
    };

} // namespace libjst

namespace seqan3 {
    template <>
    constexpr bool add_enum_bitwise_operators<libjst::node_state> = true;
} // namespace seqan3

namespace libjst {

    class node_descriptor {
    private:

        class break_descriptor {
        private:

            friend node_descriptor;

            node_state _state;

            constexpr explicit break_descriptor(node_state const state) noexcept : _state{state} {
            }
        public:

            constexpr break_descriptor() = default;

            constexpr bool from_left_begin() const noexcept {
                return node_descriptor::is_active(_state, node_state::left_begin);
            }

            constexpr bool from_left_end() const noexcept {
                return node_descriptor::is_active(_state, node_state::left_end);
            }

            constexpr bool from_right_begin() const noexcept {
                return node_descriptor::is_active(_state, node_state::right_begin);
            }

            constexpr bool from_right_end() const noexcept {
                return node_descriptor::is_active(_state, node_state::right_end);
            }
        };

        node_state _state{};
        bool _on_alternate_path{};

    public:

        constexpr node_descriptor() = default;
        constexpr explicit node_descriptor(node_state const state) noexcept {
            activate_state(state);
        }

        constexpr node_descriptor & operator=(node_state const state) noexcept {
            activate_state(state);
            return *this;
        }

        constexpr void activate_state(node_state const state) noexcept {
            _state = state;
            if (is_active(_state, node_state::variant))
                _on_alternate_path = true;
        }

        constexpr bool from_reference() const noexcept {
            return !from_variant();
        }

        constexpr bool from_variant() const noexcept {
            return is_active(_state, node_state::variant);
        }

        constexpr bool is_branching() const noexcept {
            return is_active(_state, node_state::right_begin);
        }

        constexpr bool on_alternate_path() const noexcept {
            return _on_alternate_path;
        }

        constexpr void toggle_alternate_path() noexcept {
            _on_alternate_path = true;
        }

        constexpr break_descriptor left_break() const noexcept {
            using seqan3::operator|;
            using seqan3::operator&;
            return is_left_only() ? break_descriptor{node_state::left_begin}
                                  : break_descriptor{_state & (node_state::left_begin | node_state::left_end)};
        }

        constexpr break_descriptor right_break() const noexcept {
            using seqan3::operator|;
            using seqan3::operator&;
            return is_left_only() ? break_descriptor{node_state::left_end}
                                  : break_descriptor{_state & (node_state::right_begin | node_state::right_end)};
        }

        constexpr operator node_state() const noexcept {
            return _state;
        }

        template <seqan3::cereal_output_archive archive_t>
        void save(archive_t & oarchive) const
        {
            oarchive(_state, _on_alternate_path);
        }

        template <seqan3::cereal_input_archive archive_t>
        void load(archive_t & iarchive)
        {
            iarchive(_state, _on_alternate_path);
        }

    private:

        constexpr bool is_left_only() const noexcept {
            using seqan3::operator|;
            return is_active(_state, node_state::left_begin | node_state::left_end);
        }

        static constexpr bool is_active(node_state const state, node_state const query_state) noexcept {
            using seqan3::operator&;
            return (state & query_state) == query_state;
        }

        constexpr friend bool operator==(node_descriptor const & lhs, node_descriptor const & rhs) noexcept = default;
    };

    template <typename char_t, typename char_traits_t, typename node_descriptor_t>
        requires std::same_as<std::remove_cvref_t<node_descriptor_t>, node_descriptor>
    inline std::basic_ostream<char_t, char_traits_t> & operator<<(std::basic_ostream<char_t, char_traits_t> & stream,
                                                                  node_descriptor_t && descriptor)
    {
        using namespace std::literals;

        auto to_msg = [] (node_descriptor_t const & descriptor) -> std::string_view {
            switch (static_cast<node_state>(descriptor)) {
                case node_state::branching_after_left_end: return "branching_after_left_end"sv;
                case node_state::last_branching_after_left_end: return "last_branching_after_left_end"sv;
                case node_state::branching_after_left_begin: return "branching_after_left_begin"sv;
                case node_state::last_branching_after_left_begin: return "last_branching_after_left_begin"sv;
                case node_state::non_branching_left_only: return "non_branching_left_only"sv;
                case node_state::last_non_branching_left_only: return "last_non_branching_left_only"sv;
                case node_state::non_branching_including_left: return "non_branching_including_left"sv;
                case node_state::non_branching_after_left: return "non_branching_after_left"sv;
                case node_state::variant: return "variant"sv;
                default: return "unkown"sv;
            };
        };

        stream << "<path = " << (descriptor.on_alternate_path() ? "alt"sv : "ref"sv)
                             << " state = " << to_msg(descriptor)
               << ">";
        return stream;
    }
}  // namespace libjst
