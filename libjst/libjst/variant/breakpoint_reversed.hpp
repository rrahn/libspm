// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides reverse breakpoint implementation.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/type_traits.hpp>

#include <libjst/variant/breakpoint.hpp>
#include <libjst/variant/concept.hpp>

namespace libjst
{
    class breakpoint_reversed {
    private:

        breakpoint _wrappee{};
        size_t _offset{};

    public:

        using value_type = typename breakpoint::value_type;

        constexpr breakpoint_reversed() = delete;
        ~breakpoint_reversed() = default;

        constexpr explicit breakpoint_reversed(breakpoint wrappee, size_t offset) noexcept :
            _wrappee{std::move(wrappee)},
            _offset{offset}
        {}

    private:

        template <typename me_t>
            requires std::same_as<std::remove_cvref_t<me_t>, breakpoint_reversed>
        constexpr friend value_type tag_invoke(std::tag_t<libjst::low_breakend>, me_t && me) noexcept
        {
            return me._offset - libjst::high_breakend(me._wrappee);
        }

        template <typename me_t>
            requires std::same_as<std::remove_cvref_t<me_t>, breakpoint_reversed>
        constexpr friend value_type tag_invoke(std::tag_t<libjst::high_breakend>, me_t && me) noexcept
        {
            return me._offset - libjst::low_breakend(me._wrappee);
        }

        template <typename cpo_t, typename me_t>
            requires std::same_as<std::remove_cvref_t<me_t>, breakpoint_reversed> &&
                     std::tag_invocable<cpo_t, jst::contrib::member_type_t<me_t, breakpoint>>
        constexpr friend auto tag_invoke(cpo_t cpo, me_t && me) noexcept
            -> std::tag_invoke_result_t<cpo_t, jst::contrib::member_type_t<me_t, breakpoint>>
        {
            using member_t = jst::contrib::member_type_t<me_t, breakpoint>;
            return cpo((member_t &&)me._wrappee);
        }

        constexpr friend bool operator==(breakpoint_reversed const &, breakpoint_reversed const &) noexcept =  default;
        constexpr friend std::strong_ordering operator<=>(breakpoint_reversed const &, breakpoint_reversed const &) noexcept = default;
    };

    template <typename stream_t, typename breakpoint_t>
        requires std::same_as<std::remove_cvref_t<breakpoint_t>, breakpoint_reversed>
    inline stream_t & operator<<(stream_t & stream, breakpoint_t && bp) {
        using namespace std::literals;
        return stream << "(" << libjst::low_breakend(bp) << ".." << libjst::high_breakend(bp) << "]";
    }
}  // namespace libjst
