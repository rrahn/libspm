// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides type traits utilities.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

namespace libio
{
    namespace detail
    {
        template <typename member_t, typename class_t>
        constexpr auto to_member_type(class_t const &) noexcept -> member_t class_t::*;

        template <typename class_t, typename member_t>
        using mem_fn_t = decltype(std::mem_fn(to_member_type<member_t>(std::declval<class_t &&>())));

    } // namespace detail

    template <typename class_t, typename member_t>
    using member_type_t = std::invoke_result_t<detail::mem_fn_t<class_t, member_t>, class_t &&>;

} // namespace libio
