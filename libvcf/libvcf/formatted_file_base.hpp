// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <string>
#include <seqan3/std/filesystem>

namespace libvcf {

template <typename char_t, typename char_traits_t = std::char_traits<char_t>>
struct raw_record {
    std::basic_string<char_t, char_traits_t> value{};  //!\brief The extracted raw record value.
};

template <typename format_t>
class formatted_file_base {
private:
    struct format_base;

    using char_t = char;
    using char_traits_t = std::char_traits<char_t>;
    using stream_t = std::basic_istream<char_t, char_traits_t>;
    using raw_record_t = raw_record<char_t, char_traits_t>;

    std::unique_ptr<stream_t> _stream{};
    format_t _format{};
public:

    using value_type = raw_record_t;
    // stream
    // stream needs to be opened inside of constructor

    template <typename format_t>
    formatted_file_base(std::filesystem::path && file_path) {
        // set_format()
        // format_file_path?
        // extract the extensions from the filepath
        // we do not care whether the extension of the file has the valid format.
        // valid format is checked by content of the file.

        open_stream(std::move(file_path));
        // now we need some form of type erasure for the format itself.
    }
    // handle single record and all other elements?

    value_type read_record() {
        return _format->get(*_stream);
    }

protected:

    /*!\brief Opens an input file stream with the given path.
     * \param file_path The file path.
     */
    void open_stream(std::filesystem::path && file_path) {
        _stream = std::make_unique<std::basic_ifstream<char_t, char_traits_t>>(file_path);
    }

    template <typename format_t>
    void set_format(format_t format) {
        _format = std::make_unique<format_impl<format_t>>(std::move(format));
    }

private:

    // This allows me to read into a raw_record, but we don't have the actual type for the record itself.
    struct format_base {
        virtual raw_record_t get(stream_t &);
    };

    template <typename format_t>
    struct format_impl : public format_base {
        format_t _format;

        raw_record_t get(stream_t & stream) override {
            return _format.get(stream);
        }
    };
};
} // namesapce libvcf
