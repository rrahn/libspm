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

#include <exception>
#include <tuple>

#include <libspm/closure_adaptor.hpp>
#include <libspm/execute/concept_receiver.hpp>
#include <libspm/execute/then.hpp>

namespace execute
{

    namespace _make_sender
    {

        template <typename receiver_t, typename fn_t, typename arg_tuple_t>
        struct command {
            receiver_t receiver;
            fn_t fn;
            arg_tuple_t arg_tuple;

            void start() {
                try {
                    auto value = std::apply([&](auto &&...args) {
                        return std::invoke((fn_t&&)fn, (decltype(args)&&)args...);
                    }, arg_tuple);
                    execute::set_value((receiver_t&&)receiver, std::move(value));
                } catch (...) {
                    execute::set_error((receiver_t&&)receiver, std::current_exception());
                }
            };
        };

        template <typename fn_t, typename ...args_t>
        struct sender {
            fn_t fn;
            std::tuple<args_t...> arg_tuple;

            template <typename receiver_t>
            auto connect(receiver_t&& receiver) noexcept -> command<receiver_t, fn_t, std::tuple<args_t...>> {
                using command_t = command<receiver_t, fn_t, std::tuple<args_t...>>;
                return command_t{(receiver_t&&)receiver, (fn_t&&) fn, std::move(arg_tuple)};
            }
        };

        inline struct closure
        {
            template <typename fn_t, typename ...args_t>
            auto operator()(fn_t&& fn, args_t&& ...args) const
                noexcept(std::is_nothrow_constructible_v<sender<fn_t, args_t...>>)
                -> sender<fn_t, args_t...>
            {
                return sender<fn_t, args_t...>{(fn_t&&)fn, (args_t &&)args...};
            }

            template <typename ...args_t>
            auto operator()(args_t&& ...args) const
                noexcept(noexcept(spm::make_closure(std::declval<closure>(), (args_t&&)args...)))
                -> spm::closure_result_t<closure, args_t...>
            {
                return spm::make_closure(*this, (args_t&&)args...);
            }
        } make_sender;

    } // namespace _make_sender

    using _make_sender::make_sender;

} // namespace execute
