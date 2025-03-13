// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides make stream factory.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libspm/execute/ready_done.hpp>
#include <libspm/execute/then.hpp>

#include <libspm/copyable_box.hpp>
#include <libspm/execute/concept_operation.hpp>
#include <libspm/execute/concept_receiver.hpp>
#include <libspm/execute/concept_stream.hpp>

namespace execute
{

    namespace _for_each_stream
    {

        template <typename parent_stream_t, typename fn_t>
        class sender
        {
            template <typename receiver_t>
            class command;

            template <typename receiver_t>
            class next_receiver;

            parent_stream_t _parent_stream;
            fn_t _fn;

        public:
            explicit sender(parent_stream_t parent_stream, fn_t fn) :
                _parent_stream{(parent_stream_t &&) parent_stream},
                _fn{(fn_t &&) fn}
            {
            }

            template <typename receiver_t>
            command<receiver_t> connect(receiver_t && receiver) noexcept {
                return command<receiver_t>{(parent_stream_t &&) _parent_stream,
                                           (fn_t &&)_fn,
                                           (receiver_t &&)receiver};
            }
        };

        template <typename parent_stream_t, typename fn_t>
        template <typename receiver_t>
        class sender<parent_stream_t, fn_t>::command
        {
            friend sender;

            using fn_box_t = jst::contrib::copyable_box<std::remove_reference_t<fn_t>>;

            parent_stream_t _parent_stream;
            fn_box_t _fn;
            receiver_t _receiver;
            bool _eof{false};
            std::exception_ptr _error{};


            explicit command(parent_stream_t parent_stream, fn_t fn, receiver_t receiver) noexcept :
                _parent_stream{(parent_stream_t &&) parent_stream},
                _fn{(fn_t &&) fn},
                _receiver{(receiver_t &&) receiver}
            {
            }

        public:
            void start() noexcept
            {
                try {
                    while (!eof()) {
                        auto next_sender = execute::next(_parent_stream); //| execute::then(*_fn);
                        auto next_command = execute::connect(next_sender, next_receiver<receiver_t>{*this});
                        execute::start(next_command);
                    }

                    if (_error) {
                       std::rethrow_exception(_error);
                    }

                    auto cleanup_sender = execute::cleanup(_parent_stream);
                    auto cleanup_command = execute::connect(cleanup_sender, next_receiver<receiver_t>{*this});
                    execute::start(cleanup_command);
                    execute::set_value((receiver_t &&)_receiver);
                } catch (...) {
                    auto cleanup_sender = execute::cleanup(_parent_stream);
                    auto cleanup_command = execute::connect(cleanup_sender, next_receiver<receiver_t>{*this});
                    execute::start(cleanup_command);
                    ((receiver_t &&)_receiver).set_error(std::current_exception());
                    execute::set_error((receiver_t &&)_receiver, std::current_exception());
                }
            }
        private:

            bool eof() const noexcept {
                return _eof;
            }

            void set_done() noexcept {
                _eof = true;
            }

            void set_error(std::exception_ptr error) noexcept {
                set_done();
                _error = error;
            }

        };

        template <typename parent_stream_t, typename fn_t>
        template <typename receiver_t>
        class sender<parent_stream_t, fn_t>::next_receiver
        {
            command<receiver_t> & _host;
        public:
            next_receiver(command<receiver_t> & host) noexcept : _host{host}
            {}

            template <typename ...args_t>
            void set_value(args_t&&...args) && noexcept
            {
                std::invoke(*_host._fn, (args_t&&) args...);
            }

            void set_done() && noexcept
            {
                _host.set_done();
            }

            void set_error(std::exception_ptr error) && noexcept
            {
                _host.set_error(error);
            }
        };

        inline struct closure
        {
            template <typename parent_stream_t, typename fn_t>
            auto operator()(parent_stream_t&& parent_stream, fn_t&& fn) const
                noexcept(std::is_nothrow_constructible_v<sender<parent_stream_t, fn_t>>)
                -> sender<parent_stream_t, fn_t>
            {
                return sender<parent_stream_t, fn_t>{(parent_stream_t &&) parent_stream, (fn_t &&) fn};
            }

            template <typename fn_t>
            auto operator()(fn_t && fn) const
                noexcept(noexcept(jst::contrib::make_closure(std::declval<closure>(), (fn_t&&)fn)))
                -> jst::contrib::closure_result_t<closure, fn_t>
            {
                return jst::contrib::make_closure(closure{}, (fn_t &&)fn);
            }

        } for_each_stream;

    } // namespace _for_each_stream

    using _for_each_stream::for_each_stream;
} // namespace execute
