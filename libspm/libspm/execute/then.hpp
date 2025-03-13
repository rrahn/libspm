// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides then sender.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <exception>

#include <libspm/closure_adaptor.hpp>
#include <libspm/execute/concept_receiver.hpp>
#include <libspm/execute/concept_sender.hpp>

namespace execute
{

    namespace _then
    {

        template <typename receiver_t, typename fn_t>
        struct command {
            receiver_t receiver;
            fn_t fn;

            template <typename ...args_t>
            void set_value(args_t &&...args) && {
                try {
                    if constexpr (std::same_as<std::invoke_result_t<fn_t, args_t...>, void>) {
                        std::invoke((fn_t&&)fn, (args_t&&)args...);
                        execute::set_value((receiver_t&&)receiver);
                    } else {
                        auto res = std::invoke((fn_t&&)fn, (args_t&&)args...);
                        execute::set_value((receiver_t&&)receiver, std::move(res));
                    }
                } catch (...) {
                    execute::set_error((receiver_t&&)receiver, std::current_exception());
                }
            }

            void set_done() && noexcept {
                execute::set_done((receiver_t&&)receiver);
            }

            void set_error(std::exception_ptr err) && {
                execute::set_error((receiver_t&&)receiver, std::move(err));
            }
        };

        template <typename parent_sender_t, typename fn_t>
        struct sender {
            parent_sender_t parent_sender;
            fn_t fn;

            template <typename receiver_t>
            auto connect(receiver_t && receiver) noexcept
                -> execute::operation_t<parent_sender_t, command<receiver_t, fn_t>> {
                using command_t = command<receiver_t, fn_t>;
                return execute::connect((parent_sender_t&&)parent_sender,
                                        command_t{(receiver_t&&)receiver, (fn_t&&) fn});
            }
        };

        inline struct closure
        {
            template <typename sender_t, typename fn_t>
            auto operator()(sender_t && parent_sender, fn_t && fn) const
                noexcept(std::is_nothrow_constructible_v<sender<sender_t, fn_t>>)
                -> sender<sender_t, fn_t>
            {
                return sender<sender_t, fn_t>{(sender_t &&)parent_sender, (fn_t &&) fn};
            }

            template <typename fn_t>
            auto operator()(fn_t && fn) const
                noexcept(noexcept(spm::make_closure(std::declval<closure>(), (fn_t&&)fn)))
                -> spm::closure_result_t<closure, fn_t>
            {
                return spm::make_closure(closure{}, (fn_t &&)fn);
            }
        } then;

    } // namespace _then

    using _then::then;

} // namespace execute
