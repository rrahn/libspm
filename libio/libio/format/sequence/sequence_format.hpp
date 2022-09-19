// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <iosfwd>
#include <memory>
#include <type_traits>

#include <libio/format/sequence/sequence_token.hpp>
#include <libio/utility/tag_invoke.hpp>

namespace libio
{
    template <typename ...formats_t>
    class sequence_format
    {
    private:
        struct format_base;

        template <typename format_t, typename stream_t>
        using format_token_t = tag_invoke_result_t<tag_t<libio::format_token>, format_t &, stream_t &>;

        template <typename stream_t>
        using sequence_token_t = sequence_token<format_token_t<formats_t, stream_t>...>;

        using valid_formats_t = std::vector<std::unique_ptr<format_base>>;

        static constexpr struct make_vector_helper
        {
            template <typename... base_formats_t>
            constexpr auto operator()(base_formats_t &&...formats) const
            {
                valid_formats_t vec{};
                vec.reserve(sizeof...(formats));
                emplace_next(vec, type_erase(std::forward<base_formats_t>(formats))...);
                return vec;
            }

        private:

            template <typename format_t, typename... base_formats_t>
            static constexpr auto emplace_next(valid_formats_t &vec, format_t &&format, base_formats_t &&...remaining) noexcept
            {
                vec.emplace_back(std::forward<format_t>(format));
                if constexpr (sizeof...(remaining) > 0)
                    emplace_next(vec, std::forward<base_formats_t>(remaining)...);
            }

            template <typename format_t>
            static constexpr auto type_erase(format_t &&format) noexcept
            {
                return std::make_unique<format_impl<format_t>>(std::forward<format_t>(format));
            }
        } make_unique_formats{};

        valid_formats_t _valid_formats{};
        format_base * _selected_format{};

    public:
        // we are called with a specific list of formats
        template <typename ..._formats_t>
            requires (sizeof...(_formats_t) >= 1)
        sequence_format(_formats_t && ...formats)
            : _valid_formats{make_unique_formats(std::forward<_formats_t>(formats)...)}
        {
        }

    private:

        constexpr friend bool tag_invoke(tag_t<libio::select_format>,
                                         sequence_format & me,
                                         std::filesystem::path const & path) noexcept
        {
            me._selected_format = nullptr;
            std::ranges::for_each(me._valid_formats, [&] (auto const & format_ptr)
            {
                auto const & extensions = format_ptr->valid_extensions();
                if (auto it = std::ranges::find(extensions, static_cast<std::string>(path.extension())); it != std::ranges::end(extensions))
                    me._selected_format = format_ptr.get();
            });
            return me._selected_format != nullptr;
        }

        template <typename char_t, typename char_traits_t>
        constexpr friend sequence_token_t<std::basic_istream<char_t, char_traits_t>>
        tag_invoke(tag_t<libio::format_token>,
                   sequence_format const & fmt,
                   std::basic_istream<char_t, char_traits_t> & istream)
        {
            assert(fmt._selected_format != nullptr);
            return fmt._selected_format->format_token(istream);
        }


        using common_extensions_t =
            std::common_reference_t<decltype(libio::valid_extensions(std::declval<formats_t const &>()))...>;
        // type erased format
        struct format_base
        {
            virtual ~format_base() = default;
            virtual sequence_token_t<std::istream> format_token(std::istream &) const = 0;
            virtual common_extensions_t valid_extensions() const = 0;
        };

        template <typename format_t>
        struct format_impl final : public format_base
        {
            format_t _format{}; // the specific format for which we get the token?

            format_impl(format_t format) : _format{std::forward<format_t>(format)}
            {
            }
            virtual ~format_impl() = default;

            sequence_token_t<std::istream> format_token(std::istream &istream) const override
            {
                return libio::format_token(_format, istream);
            }

            common_extensions_t valid_extensions() const override
            {
                return libio::valid_extensions(_format);
            }
        };
    };

    template <typename ...formats_t>
    sequence_format(formats_t &&...) -> sequence_format<formats_t...>;
} // namespace libio
