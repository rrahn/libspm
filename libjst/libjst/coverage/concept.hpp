// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concepts for breakpoint coverages.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/std/tag_invoke.hpp>

namespace libjst
{

    namespace _coverage_intersection {
        inline constexpr struct _cpo  {
            template <typename coverage1_t, typename coverage2_t>
                requires std::tag_invocable<_cpo, coverage1_t, coverage2_t>
            constexpr auto operator()(coverage1_t && c1, coverage2_t && c2) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, coverage1_t, coverage2_t>)
                -> std::tag_invoke_result_t<_cpo, coverage1_t, coverage2_t>
            {
                return std::tag_invoke(_cpo{}, (coverage1_t &&) c1, (coverage2_t &&) c2);
            }
        } coverage_intersection;
    } // namespace _coverage_intersection
    using _coverage_intersection::coverage_intersection;

    namespace _coverage_difference {
        inline constexpr struct _cpo  {
            template <typename coverage1_t, typename coverage2_t>
                requires std::tag_invocable<_cpo, coverage1_t, coverage2_t>
            constexpr auto operator()(coverage1_t && c1, coverage2_t && c2) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, coverage1_t, coverage2_t>)
                -> std::tag_invoke_result_t<_cpo, coverage1_t, coverage2_t>
            {
                return std::tag_invoke(_cpo{}, (coverage1_t &&) c1, (coverage2_t &&) c2);
            }
        } coverage_difference;
    } // namespace _coverage_difference
    using _coverage_difference::coverage_difference;

    namespace _get_domain {
        inline constexpr struct _cpo  {
            template <typename coverage_t>
                requires std::tag_invocable<_cpo, coverage_t>
            constexpr auto operator()(coverage_t && coverage) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, coverage_t>)
                -> std::tag_invoke_result_t<_cpo, coverage_t>
            {
                return std::tag_invoke(_cpo{}, (coverage_t &&) coverage);
            }

        private:

            template <typename coverage_t>
                requires requires (coverage_t coverage) { { std::forward<coverage_t &&>(coverage).get_domain() }; }
            constexpr friend auto tag_invoke(_cpo, coverage_t && coverage)
                noexcept(noexcept(std::declval<coverage_t &&>().get_domain()))
                -> decltype(std::declval<coverage_t &&>().get_domain())
            {
                return ((coverage_t &&)coverage).get_domain();
            }
        } get_domain;
    } // namespace _get_domain
    using _get_domain::get_domain;

    template <typename coverage_t>
    using coverage_domain_t = std::remove_cvref_t<std::tag_invoke_result_t<_get_domain::_cpo, coverage_t>>;
}  // namespace libjst
