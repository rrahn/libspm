// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libio::token_get_area.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <ranges>

#include <libio/file/tokenizer_streambuffer_adaptor.hpp>
#include <libio/file/consume_tokenizer.hpp>
#include <libio/file/until_tokenizer.hpp>

namespace libio
{
    // wrapper around
    template <typename tokenizer_t>
    class token_get_area : public tokenizer_t // until_tokenizer<tokenizer_streambuffer_adaptor<streambuffer_t>>
    {
        // using token_buffer_t = tokenizer_streambuffer_adaptor<streambuffer_t>;
        using base_t = tokenizer_t;
    public:

        using typename base_t::char_type;
        using typename base_t::traits_type;
        using typename base_t::int_type;
        using typename base_t::pos_type;
        using typename base_t::off_type;

        token_get_area() = default;
        token_get_area(token_get_area const &) = delete;
        token_get_area(token_get_area &&) = default;
        token_get_area &operator=(token_get_area const &) = delete;
        token_get_area &operator=(token_get_area &&) = default;

        template <typename stream_buffer_t, typename token_t>
            requires std::invocable<token_t &&, tokenizer_streambuffer_adaptor<stream_buffer_t>>
        constexpr token_get_area(stream_buffer_t * stream_buffer, token_t && token) :
            tokenizer_t{std::invoke((token_t &&)token, tokenizer_streambuffer_adaptor{stream_buffer})}
        {
        }

        ~token_get_area() = default;

        constexpr std::ranges::iterator_t<base_t> begin() noexcept
        {
            return base_t::begin();
        }
        void begin() const noexcept = delete;

        constexpr std::ranges::sentinel_t<base_t> end() noexcept
        {
            return base_t::end();
        }

        void end() const noexcept = delete;
    };

    template <typename streambuffer_t, typename token_t>
        requires std::invocable<token_t &&, tokenizer_streambuffer_adaptor<streambuffer_t>>
    token_get_area(streambuffer_t *, token_t &&)
        -> token_get_area<std::invoke_result_t<token_t &&, tokenizer_streambuffer_adaptor<streambuffer_t>>>;
}  // namespace libio
