// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides utilities for fast tokenization.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <ranges>

#include <libio/utility/tag_invoke.hpp>

namespace libio
{
    namespace _read_token {

        inline constexpr struct _cpo {

            template <typename value_t, typename irange_t>
                requires requires (value_t & value) { std::back_inserter(value); }
            friend void tag_invoke(_cpo, value_t & value, irange_t & chunked_buffer)
            {
                for (auto && chunk : chunked_buffer)
                {
                    value.reserve(value.capacity() + std::ranges::size(chunk));
                    std::ranges::copy(chunk, std::back_inserter(value));
                }
            }

            template <typename value_t, typename irange_t>
                requires std::ranges::input_range<irange_t>
            auto operator()(value_t & value, irange_t & irange) const
                noexcept(is_nothrow_tag_invocable_v<_cpo, value_t &, irange_t &>)
                -> tag_invoke_result_t<_cpo, value_t &, irange_t &>
            {
                return libio::tag_invoke(_cpo{}, value, irange);
            }
        } read_token{};
    } // namespace _read_token

    using _read_token::read_token;


    namespace _detokenize {
        inline constexpr struct _cpo {

            template <typename token_t, typename value_t>
            friend void tag_invoke(_cpo,
                                   [[maybe_unused]] token_t const & token,
                                   [[maybe_unused]] value_t const & value)
            { // in default action we say there is no overload found!
                // maybe we can set error on token?
                // set_error(token, error_code{...}); // then the token can decide wether to throw or to set an ec.
                throw std::runtime_error{"No known overload found for detokenize."};
            }

            template <typename token_t, typename value_t>
            auto operator()(token_t && token, value_t && value) const
                noexcept(is_nothrow_tag_invocable_v<_cpo, token_t &&, value_t &&>)
                -> tag_invoke_result_t<_cpo, token_t &&, value_t &&>
            {
                return tag_invoke(_cpo{}, (token_t &&)token, (value_t&&)value);
            }
        } detokenize_to;
    } // namespace _detokenize
    using _detokenize::detokenize_to;

}  // namespace libio
