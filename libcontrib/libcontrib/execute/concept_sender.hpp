// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides sender concepts.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/std/tag_invoke.hpp>

namespace execute
{
    namespace _connect {
        inline constexpr struct _cpo  {
            template <typename sender_t, typename receiver_t>
                requires std::tag_invocable<_cpo, sender_t, receiver_t>
            constexpr auto operator()(sender_t && sender, receiver_t && receiver) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, sender_t, receiver_t>)
                -> std::tag_invoke_result_t<_cpo, sender_t, receiver_t>
            {
                return std::tag_invoke(_cpo{}, (sender_t &&)sender, (receiver_t &&)receiver);
            }

        private:

            template <typename sender_t, typename receiver_t>
                requires requires (sender_t sender, receiver_t receiver) {
                    { std::forward<sender_t &&>(sender).connect((receiver_t &&)receiver) };
                }
            constexpr friend auto tag_invoke(_cpo, sender_t && sender, receiver_t && receiver)
                noexcept(noexcept(std::declval<sender_t &&>().connect((receiver_t &&)receiver)))
                -> decltype(std::declval<sender_t &&>().connect((receiver_t &&)receiver))
            {
                return ((sender_t &&)sender).connect((receiver_t &&)receiver);
            }
        } connect;
    } // namespace _connect
    using _connect::connect;

    template <typename sender_t, typename receiver_t>
    using operation_t = std::tag_invoke_result_t<std::tag_t<execute::connect>, sender_t, receiver_t>;
}  // namespace execute
