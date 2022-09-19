// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concept defintion and CPOs for the formats.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <filesystem>
#include <ranges>
#include <string>
#include <vector>

#include <libio/utility/tag_invoke.hpp>

namespace libio
{
    namespace _valid_extensions {

        inline constexpr struct _cpo {

            template <typename format_t>
                requires requires (format_t const & f) { { f.valid_extensions() } -> std::same_as<std::vector<std::string> const &>; }
            friend auto tag_invoke(_cpo, format_t const & format)
                noexcept(noexcept(format.valid_extensions()))
                -> decltype(std::declval<format_t const &>().valid_extensions())
            {
                return format.valid_extensions();
            }

            // Delegate to the friend overload if default is not found
            template <typename format_t>
                requires tag_invocable<_cpo, format_t>
            auto operator()(format_t && format) const
                noexcept(is_nothrow_tag_invocable_v<_cpo, format_t &&>)
                -> tag_invoke_result_t<_cpo, format_t &&>
            {
                return libio::tag_invoke(_cpo{}, std::forward<format_t>(format));
            }
        } valid_extensions{};

    } // namespace _valid_extensions

    using _valid_extensions::valid_extensions; //!\brief CPO defintion to get the format extensions.

    namespace _select_format {

        inline constexpr struct _cpo {

            template <typename format_t>
                requires requires (format_t const & f) { { libio::valid_extensions(f) } -> std::ranges::input_range; }
            friend bool tag_invoke(_cpo, format_t const & format, std::filesystem::path const & path)
            {
                auto const & extensions = libio::valid_extensions(format);
                return std::ranges::find_if(extensions, path.extension()) != std::ranges::end(extensions);
            }

            // Delegate to the friend overload if default is not found
            template <typename format_t>
                requires tag_invocable<_cpo, format_t, std::filesystem::path const &>
            auto operator()(format_t && format, std::filesystem::path const & path) const
                noexcept(is_nothrow_tag_invocable_v<_cpo, format_t &&, std::filesystem::path const &>)
                -> tag_invoke_result_t<_cpo, format_t &&, std::filesystem::path const &>
            {
                return libio::tag_invoke(_cpo{}, (format_t &&)format, path);
            }
        } select_format{};

    } // namespace _select_format

    using _select_format::select_format; //!\brief CPO defintion to get an unformatted record.

    namespace _format_token {

        inline constexpr struct _cpo {
            // TODO: Is there a meaningful default?

            // Delegate to the friend overload if default is not found
            template <typename format_t, typename stream_t>
                requires tag_invocable<_cpo, format_t, stream_t &>
            auto operator()(format_t && format, stream_t & stream) const
                noexcept(is_nothrow_tag_invocable_v<_cpo, format_t &&, stream_t &>)
                -> tag_invoke_result_t<_cpo, format_t &&, stream_t &>
            {
                return libio::tag_invoke(_cpo{}, std::forward<format_t>(format), stream);
            }
        } format_token{};

    } // namespace _format_token

    using _format_token::format_token; //!\brief CPO defintion to read a record.

    // When is this a valid input format?
    // We need to be able to call read record on it, together with a stream.
    // template <typename format_t, typename stream_t>
    // concept file_format_input = requires (format_t & f, stream_t & s) {
    //     { read_record(f, s) }; // What is the return type?
    // };
}  // namespace libio
