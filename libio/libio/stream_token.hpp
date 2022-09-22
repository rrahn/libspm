// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of the stream token.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <iosfwd>

#include <libio/file/consume_tokenizer.hpp>
#include <libio/file/token_get_area.hpp>
#include <libio/file/tokenization.hpp>
#include <libio/utility/tag_invoke.hpp>

namespace libio
{
    template <typename stream_t, typename token_t>
    class stream_token
    {
    private:
        using char_t = typename stream_t::char_type;
        using traits_t = typename stream_t::traits_type;
        // we already defined the range?
        using buffer_t = decltype(token_get_area{std::declval<stream_t &>().rdbuf(), std::declval<token_t &&>()});
        using get_area_t = consume_tokenizer<buffer_t>;
    public:

        // when we create, we want to move to the next record already
        stream_t *_stream{};
        get_area_t _get_area;

        // not sure how we handle this efficiently!
        stream_token(stream_t &stream, token_t &&token) :
            _stream{std::addressof(stream)},
            _get_area{std::in_place_type<buffer_t>, _stream->rdbuf(), (token_t &&)token}
        {
        }
        stream_token(stream_token const &) = delete;
        stream_token(stream_token && other) noexcept :
            _stream{other._stream},
            _get_area{std::move(other._get_area)}
        {
            other._stream = nullptr;
        }
        stream_token &operator=(stream_token const &) = delete;
        stream_token &operator=(stream_token && other) noexcept
        {
            _stream = other._stream;
            other._stream = nullptr;
            _get_area = std::move(other._get_area);
            return *this;
        }
        ~stream_token()
        {
            if (_stream != nullptr)
            {
                // ensure that the get area is consumed for this token!
                for (auto it = _get_area.begin(); it != _get_area.end(); ++it)
                {}
                // now we need to check if the stream is also empty!
                if (_stream->rdbuf()->sgetc() == traits_t::eof())
                    _stream->setstate(std::ios_base::eofbit);
            }
        }

    protected:
        get_area_t &get_area() noexcept
        {
            return _get_area;
        }

        get_area_t const &get_area() const noexcept
        {
            return _get_area;
        }

    private:
        // what does this mean? how can materialise this now?
        // why not having the record handling the tokenization/detokenization
        // template <typename record_t>
        // friend auto tag_invoke(tag_t<libio::tokenize>, stream_token &me, record_t &record)
        // {
        // }

        // in this form we can handle the stream and run the options to get the stream.
        // not used as default?
        // template <typename value_t>
        // friend auto tag_invoke(tag_t<libio::detokenize_to>, stream_token &token, value_t &value)
        // {
        //     libio::read_token(value, token.get_area());
        // }
    };

    template <typename stream_t, typename token_t>
    stream_token(stream_t &, token_t) -> stream_token<stream_t, token_t>;
} // namespace libio
