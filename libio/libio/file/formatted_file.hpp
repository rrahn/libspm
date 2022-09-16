// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <string>
#include <seqan3/std/filesystem>

#include <seqan/stream.h>

#include <libio/file/tokenized_stream.hpp>
#include <libio/format/format_concept.hpp>
#include <libio/utility/tag_invoke.hpp>


namespace libio
{
    template <typename format_t>
    class formatted_file
    {
    private:

        // Reusing seqan base implementation for streams etc.
        using stream_t = seqan::FormattedFile<format_tag_t<format_t>, seqan::Input>;

        // where is the token coming from?
        using token_t = typename format_t::token_type;
        using tokenized_stream_t = tokenized_stream<token_t>

        // format_t:
            // exposes type of header
            // exposes type of record
            // exposes valid file extensions
            // what


        format_t _format{};
        std::unique_ptr<stream_t> _stream{};

        // every file has its own formatation rule
        // this also means we can set at every type another record type?
        // how can we orchestrate this via the format implementation?

    public:
        formatted_file(std::filesystem::path &&file_path, format_t format = {}) : _format{std::move(format)}
        {
            // Step 1: open stream (including decompression/compression)
            open_stream(std::move(file_path));
            // Step 2: use the opened stream buffer to initialise the tokenized_stream

            // Step 3:
            // format the tokens on access
        }

        // we can always return a header, or keep it disabled.
        // auto header() {
        //     libio::parse_header()
        // }

        // TODO: Specify the API!
        auto read_record() noexcept(is_nothrow_tag_invocable_v<tag_t<libio::format_record>, format_t &, stream_t &>)
            -> tag_invoke_result_t<tag_t<libio::format_record>, format_t &, stream_t &>
        {
            assert(_stream != nullptr);

            // Check applicability of a header.
            // now what are we getting?
            return libio::format_record(_format, _stream->iter);
        }

    protected:
        /*!\brief Opens an input file stream with the given path.
         * \param file_path The file path.
         */
        void open_stream(std::filesystem::path &&file_path)
        {
            _stream = std::make_unique<stream_t>(file_path.c_str());
        }

    private:
    };

    template <typename format_t>
    formatted_file(std::filesystem::path &&, format_t) -> formatted_file<format_t>;
} // namesapce libio
