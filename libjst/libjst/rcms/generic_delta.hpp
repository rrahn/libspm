// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a generic sequence variant encoding as input for the rcms object.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>

#include <libcontrib/type_traits.hpp>
#include <libcontrib/std/tag_invoke.hpp>

#include <libjst/variant/breakpoint.hpp>

namespace libjst
{

    template <typename source_t, typename coverage_t>
    class generic_delta {

        breakpoint _breakpoint{};
        source_t _alt_sequence{};
        coverage_t _coverage{};
    public:
        constexpr generic_delta() = default;
        explicit constexpr generic_delta(breakpoint breakpoint, source_t alt_sequence, coverage_t coverage) :
            _breakpoint{std::move(breakpoint)},
            _alt_sequence{std::move(alt_sequence)},
            _coverage{std::move(coverage)}
        {}

        // ----------------------------------------------------------------------------
        // Serialisation
        // ----------------------------------------------------------------------------

        template <typename archive_t>
        void load(archive_t & iarchive)
        {
            iarchive(_breakpoint, _alt_sequence, _coverage);
        }

        template <typename archive_t>
        void save(archive_t & oarchive) const
        {
            oarchive(_breakpoint, _alt_sequence, _coverage);
        }


    private:

        template <typename cpo_t, typename me_t>
            requires std::same_as<std::remove_cvref_t<me_t>, generic_delta> &&
                     std::tag_invocable<cpo_t, jst::contrib::member_type_t<me_t, breakpoint>>
        friend constexpr auto tag_invoke(cpo_t cpo, me_t && me)
            noexcept(std::is_nothrow_tag_invocable_v<cpo_t, jst::contrib::member_type_t<me_t, breakpoint>>)
            -> std::tag_invoke_result_t<cpo_t, jst::contrib::member_type_t<me_t, breakpoint>>
        {
            using fwd_breakpoint_t = jst::contrib::member_type_t<me_t, breakpoint>;
            return std::tag_invoke(cpo, static_cast<fwd_breakpoint_t>(me._breakpoint));
        }

        template <typename me_t>
            requires std::same_as<std::remove_cvref_t<me_t>, generic_delta>
        friend constexpr auto tag_invoke(std::tag_t<libjst::get_breakpoint>, me_t && me) noexcept
            -> jst::contrib::member_type_t<me_t, breakpoint>
        {
            using fwd_breakpoint_t = jst::contrib::member_type_t<me_t, breakpoint>;
            return static_cast<fwd_breakpoint_t>(me._breakpoint);
        }

        template <typename me_t>
            requires std::same_as<std::remove_cvref_t<me_t>, generic_delta>
        friend constexpr auto tag_invoke(std::tag_t<libjst::coverage>, me_t && me) noexcept
            -> jst::contrib::member_type_t<me_t, coverage_t>
        {
            using fwd_coverage_t = jst::contrib::member_type_t<me_t, coverage_t>;
            return static_cast<fwd_coverage_t>(me._coverage);
        }

        template <typename me_t>
            requires std::same_as<std::remove_cvref_t<me_t>, generic_delta>
        friend constexpr auto tag_invoke(std::tag_t<libjst::alt_sequence>, me_t && me) noexcept
            -> jst::contrib::member_type_t<me_t, source_t>
        {
            using fwd_source_t = jst::contrib::member_type_t<me_t, source_t>;
            return static_cast<fwd_source_t>(me._alt_sequence);
        }
    };
}  // namespace libjst
