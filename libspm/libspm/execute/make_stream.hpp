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

#include <libspm/execute/concept_receiver.hpp>
#include <libspm/execute/ready_done.hpp>

namespace execute
{

    namespace _make_stream
    {

        template <typename range_t>
        class stream
        {

            using iterator_t = std::ranges::iterator_t<range_t>;
            using sentinel_t = std::ranges::sentinel_t<range_t>;

            class next_sender;
            template <typename receiver_t>
            class command;

            range_t _range;
            iterator_t _current{};
            sentinel_t _end{};

        public:
            explicit stream(range_t range) : _range{(range_t &&) range}
            {
                _current = std::ranges::begin(_range);
                _end = std::ranges::end(_range);
            }

            next_sender next() noexcept
            {
                return next_sender{*this};
            }

            ready_done_sender cleanup() noexcept
            {
                return {};
            }
        };

        template <typename range_t>
        class stream<range_t>::next_sender
        {
        private:
            friend stream;

            stream &_stream;

            explicit next_sender(stream &stream) noexcept : _stream{stream}
            {
            }

        public:
            template <typename receiver_t>
            command<receiver_t> connect(receiver_t &&receiver) noexcept
            {
                return command<receiver_t>{_stream, (receiver_t &&) receiver};
            }

            // template <typename receiver_t>
            // command<receiver_t> connect(receiver_t &&) const & noexcept = delete;
        };

        template <typename range_t>
        template <typename receiver_t>
        class stream<range_t>::command
        {

            friend stream;

            stream &_stream;
            receiver_t _receiver;

            explicit command(stream &stream, receiver_t receiver) noexcept : _stream{stream},
                                                                             _receiver{(receiver_t &&) receiver}
            {
            }

        public:
            void start() noexcept
            {
                if (_stream._current == _stream._end) {
                    execute::set_done((receiver_t &&) _receiver);
                } else {
                    execute::set_value((receiver_t &&) _receiver, *_stream._current);
                    ++_stream._current;
                }
            }
        };

        inline struct closure
        {
            template <typename range_t>
            auto operator()(range_t &&range) noexcept(std::is_nothrow_constructible_v<stream<range_t>>)
                -> stream<range_t>
            {
                return stream<range_t>{(range_t &&) range};
            }

        } make_stream;

    } // namespace _make_stream

    using _make_stream::make_stream;
} // namespace execute
