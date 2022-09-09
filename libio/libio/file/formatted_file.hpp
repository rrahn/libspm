// // -----------------------------------------------------------------------------------------------------
// // Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// // Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// // This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// // shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// // -----------------------------------------------------------------------------------------------------

// #pragma once

// #include <string>
// #include <seqan3/std/filesystem>

// #include <seqan/stream.h>

// #include <libio/format/format_concept.hpp>
// #include <libio/utility/tag_invoke.hpp>

// namespace libio
// {
//     template <typename format_t>
//     class raw_file_base
//     {
//     private:
//         // Reusing seqan base implementation for streams etc.
//         using stream_t = seqan::FormattedFile<format_tag_t<format_t>, seqan::Input>;

//         format_t _format{};
//         std::unique_ptr<stream_t> _stream{};

//     public:
//         raw_file_base(std::filesystem::path &&file_path, format_t format = {}) : _format{std::move(format)}
//         {
//             // Now we just set the stream here!
//             open_stream(std::move(file_path));
//         }

//         // TODO: Specify the API!
//         auto read_record() noexcept(is_nothrow_tag_invocable_v<tag_t<libio::read_record>, format_t &, stream_t &>)
//             -> tag_invoke_result_t<tag_t<libio::read_record>, format_t &, stream_t &>
//         {
//             assert(_stream != nullptr);

//             // Check applicability of a header.
//             return libio::read_record(_format, _stream->iter);
//         }

//     protected:
//         /*!\brief Opens an input file stream with the given path.
//          * \param file_path The file path.
//          */
//         void open_stream(std::filesystem::path &&file_path)
//         {
//             _stream = std::make_unique<stream_t>(file_path.c_str());
//         }

//     private:
//     };

//     template <typename format_t>
//     raw_file_base(std::filesystem::path &&, format_t) -> raw_file_base<format_t>;
// } // namesapce libio
