// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implements the line buffer.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/utility/char_operations/predicate.hpp>

#include <libio/file/until_tokenizer.hpp>

namespace libio
{
    template <typename tokenizer_t>
    class line_tokenizer : public until_tokenizer<tokenizer_t>
    {
    private:
        using base_t = until_tokenizer<tokenizer_t>;

        // make more versatile by handling it as a multicharacter matcher.
        static constexpr auto is_newline = seqan3::is_char<'\n'> || seqan3::is_char<'\r'>;
    public:
        using typename base_t::char_type;
        using typename base_t::traits_type;
        using typename base_t::int_type;
        using typename base_t::pos_type;
        using typename base_t::off_type;

        line_tokenizer() = delete;
        constexpr explicit line_tokenizer(tokenizer_t tokenizer) noexcept : base_t{(tokenizer_t &&)tokenizer, is_newline}
        {
        }
        line_tokenizer(line_tokenizer const &) = delete;
        line_tokenizer(line_tokenizer &&) = default;

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

    template <typename tokenizer_t>
    line_tokenizer(tokenizer_t &&) -> line_tokenizer<tokenizer_t>;
} // namespace libio
