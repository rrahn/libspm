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

#include <libio/file/token_get_area.hpp>
#include <libio/file/tokenization.hpp>
#include <libio/utility/tag_invoke.hpp>

namespace libio
{

    template <typename streambuffer_t>
    class stream_token
    {
    private:
        // we already defined the range?
        using get_area_t = token_get_area<streambuffer_t>;
    public:

        // when we create, we want to move to the next record already
        get_area_t _get_area{};

        // not sure how we handle this efficiently!
        template <typename delimiter_t>
        explicit stream_token(streambuffer_t *stream_buffer, delimiter_t delim) :
            _get_area{stream_buffer, std::move(delim)}
        {
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
        template <typename value_t>
        friend auto tag_invoke(tag_t<libio::detokenize_to>, stream_token &token, value_t &value)
        {
            libio::read_token(value, token.get_area());
        }
    };

    template <typename stream_buffer_t, typename delimiter_t>
    stream_token(stream_buffer_t &, delimiter_t &&) -> stream_token<stream_buffer_t>;
} // namespace libio
