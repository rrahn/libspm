// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of the formatted stream.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iostream>

#include <libio/format/format_concept.hpp>
#include <libio/file/tokenization.hpp>

namespace libio
{
    template <typename format_t, typename stream_t = std::istream>
    class formatted_stream final : protected stream_t
    {
    private:

        format_t *_format{};
    public:
        formatted_stream() = default;
        explicit formatted_stream(format_t &format) : stream_t{}, _format{std::addressof(format)}
        {}
        explicit formatted_stream(format_t &format, stream_t && stream) :
            stream_t{std::move(stream)},
            _format{std::addressof(format)}
        {
            auto fmt_meta_tkn = libio::get_meta_token(*_format, static_cast<stream_t &>(*this));
            libio::detokenize_to(fmt_meta_tkn, *_format); // write tokenized into format.
        }
        formatted_stream(formatted_stream const &) = delete;
        formatted_stream(formatted_stream &&) = default; // does not work.

        using stream_t::rdbuf;
        using stream_t::rdstate;
        using stream_t::setstate;
        using stream_t::eof;

        auto get() -> decltype(libio::format_token(std::declval<format_t &>(), std::declval<stream_t &>()))
        {
            if (_format == nullptr || rdbuf() == nullptr)
            {
                setstate(std::ios_base::failbit);
            }
            return libio::format_token(*_format, static_cast<stream_t &>(*this));
        }
    };

    template <typename format_t>
    formatted_stream(format_t const &) -> formatted_stream<format_t>;

    template <typename format_t, typename stream_t>
    formatted_stream(format_t const &, stream_t &&) -> formatted_stream<format_t, stream_t>;

    template <typename format_t, typename stream_t, typename record_t>
        // requires format_t is detokenizable to record.
    formatted_stream<format_t, stream_t> &
    operator>> (formatted_stream<format_t, stream_t> & formatted_istream, record_t & record)
    {
        auto tkn = formatted_istream.get();
        libio::detokenize_to(tkn, record);
        return formatted_istream;
    }

} // namespace libio
