// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides operation concept.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/std/tag_invoke.hpp>

namespace execute
{
    namespace _start {
        inline constexpr struct _cpo  {
            template <typename operation_t>
                requires std::tag_invocable<_cpo, operation_t>
            constexpr auto operator()(operation_t && operation) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, operation_t>)
                -> std::tag_invoke_result_t<_cpo, operation_t>
            {
                return std::tag_invoke(_cpo{}, (operation_t &&) operation);
            }

        private:

            template <typename operation_t>
                requires requires (operation_t operation) { { std::forward<operation_t &&>(operation).start() }; }
            constexpr friend auto tag_invoke(_cpo, operation_t && operation)
                noexcept(noexcept(std::declval<operation_t &&>().start()))
                -> decltype(std::declval<operation_t &&>().start())
            {
                return ((operation_t &&)operation).start();
            }
        } start;
    } // namespace _start
    using _start::start;
}  // namespace execute
