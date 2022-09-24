// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <string>
#include <fstream>
#include <seqan3/std/filesystem>

#include <libio/file/formatted_stream.hpp>
#include <libio/format/format_concept.hpp>
#include <libio/utility/tag_invoke.hpp>


namespace libio
{
    template <typename record_t, typename format_t>
        // check compatibility with record.
    class formatted_file
    {
    private:

        class iterator;
        using sentinel = std::default_sentinel_t;

        using stream_t = std::ifstream;
        using formatted_stream_t = formatted_stream<format_t, stream_t>;


        format_t _format{};
        std::unique_ptr<formatted_stream_t> _stream{};
        record_t _cached_record{};
        bool _is_eof{false};

        // every file has its own formatation rule
        // this also means we can set at every type another record type?
        // how can we orchestrate this via the format implementation?

    public:
        formatted_file(std::filesystem::path const &file_path, format_t format = {}) : _format{std::move(format)}
        {
            // Step 1: select format based on file
            libio::select_format(_format, file_path);
            // Step 2: open stream (including decompression/compression)
            open_stream(file_path);
        }

        format_t const & format() const noexcept
        {
            return _format;
        }

        // we can always return a header, or keep it disabled.
        // auto header() {
        //     libio::parse_header()
        // }

        iterator begin() noexcept
        {
            return iterator{this};
        }

        sentinel end() noexcept
        {
            return {};
        }

    protected:
        /*!\brief Opens an input file stream with the given path.
         * \param file_path The file path.
         */
        void open_stream(std::filesystem::path const &file_path)
        {
            _stream = std::make_unique<formatted_stream_t>(_format, stream_t{file_path});
        }
    };

    // template <typename format_t>
    // formatted_file(std::filesystem::path &&, format_t) -> formatted_file<format_t>;

    template <typename record_t, typename format_t>
        // check compatibility with record.
    class formatted_file<record_t, format_t>::iterator
    {
    private:

        formatted_file * _host{};

    public:
        using value_type = record_t; // return a new token abstraction
        using reference = record_t const &;
        using difference_type = std::ptrdiff_t;
        using pointer = std::add_pointer_t<record_t>;
        using category = std::input_iterator_tag;

        iterator() = default;
        constexpr explicit iterator(formatted_file *host) noexcept : _host{host}
        {
            ++(*this);  // what if empty stream?
        }

        constexpr reference operator*() const noexcept
        {
            return _host->_cached_record;
        }

        constexpr pointer operator->() const noexcept
        {
            return std::addressof(_host->_cached_record);
        }

        constexpr iterator &operator++()
        {
            _host->_is_eof = _host->_stream->eof();
            if (!_host->_is_eof) {
                _host->_cached_record.clear(); //= record_t{}; // clear record;
                *(_host->_stream) >> _host->_cached_record;
            }
            return *this;
        }

        constexpr void operator++(int)
        {
            ++(*this);
        }

        constexpr bool operator==(sentinel const &) const noexcept
        {
            return _host->_is_eof;
        }
    };
} // namesapce libio
