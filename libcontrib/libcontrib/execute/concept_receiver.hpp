// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concepts for the receiver.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/std/tag_invoke.hpp>

namespace execute
{

    namespace _set_value {
        inline constexpr struct _cpo  {
            template <typename receiver_t, typename ...values_t>
                requires std::tag_invocable<_cpo, receiver_t, values_t...>
            constexpr void operator()(receiver_t && receiver, values_t &&... values) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, receiver_t, values_t...>)
            {
                return std::tag_invoke(_cpo{}, (receiver_t &&) receiver, (values_t &&)values...);
            }

        private:

            template <typename receiver_t, typename ...values_t>
                requires requires (receiver_t receiver, values_t &&... values) {
                    { std::forward<receiver_t &&>(receiver).set_value((values_t &&)values...) };
                }
            constexpr friend void tag_invoke(_cpo, receiver_t && receiver, values_t &&... values)
                noexcept(noexcept(std::declval<receiver_t &&>().set_value((values_t &&)values...)))
            {
                return ((receiver_t &&)receiver).set_value((values_t &&)values...);
            }
        } set_value;
    } // namespace _set_value
    using _set_value::set_value;

    namespace _set_done {
        inline constexpr struct _cpo  {
            template <typename receiver_t>
                requires std::tag_invocable<_cpo, receiver_t>
            constexpr void operator()(receiver_t && receiver) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, receiver_t>)
            {
                return std::tag_invoke(_cpo{}, (receiver_t &&) receiver);
            }

        private:

            template <typename receiver_t>
                requires requires (receiver_t receiver) { { std::forward<receiver_t &&>(receiver).set_done() }; }
            constexpr friend void tag_invoke(_cpo, receiver_t && receiver)
                noexcept(noexcept(std::declval<receiver_t &&>().set_done()))
            {
                return ((receiver_t &&)receiver).set_done();
            }
        } set_done;
    } // namespace _set_done
    using _set_done::set_done;

    namespace _set_error {
        inline constexpr struct _cpo  {
            template <typename receiver_t, typename error_t>
                requires std::tag_invocable<_cpo, receiver_t, error_t>
            constexpr void operator()(receiver_t && receiver, error_t && error) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, receiver_t, error_t>)
            {
                return std::tag_invoke(_cpo{}, (receiver_t &&) receiver, (error_t &&) error);
            }

        private:

            template <typename receiver_t, typename error_t>
                requires requires (receiver_t receiver, error_t error) {
                    { std::forward<receiver_t &&>(receiver).set_error((error_t &&)error) };
                }
            constexpr friend void tag_invoke(_cpo, receiver_t && receiver, error_t && error)
                noexcept(noexcept(std::declval<receiver_t &&>().set_error((error_t &&) error)))
            {
                return ((receiver_t &&)receiver).set_error((error_t &&) error);
            }
        } set_error;
    } // namespace _set_error
    using _set_error::set_error;
}  // namespace execute
