// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides CPOs for the records.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libio/utility/tag_invoke.hpp>

namespace libio
{
    namespace _set_field
    {
        inline constexpr struct _cpo {

            template <typename record_t, typename field_code_t, typename ibuffer_t>
                requires tag_invocable<_cpo, record_t &, field_code_t const &, ibuffer_t &&>
            constexpr auto operator()(record_t & record, field_code_t const & field_code, ibuffer_t && ibuffer) const
                noexcept(is_nothrow_tag_invocable_v<_cpo, record_t &, field_code_t const &, ibuffer_t &&>)
                -> tag_invoke_result_t<_cpo, record_t &, field_code_t const &, ibuffer_t &&>
            {
                return libio::tag_invoke(_cpo{}, record, field_code, std::forward<ibuffer_t>(ibuffer));
            }

        } set_field{};
    } // namespace _set_field

    using _set_field::set_field;

}  // namespace libio
