// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implements the token for the fastq format.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <memory>
#include <variant>

#include <libio/file/tokenization.hpp>

namespace libio
{
    template <typename ...tokens_t>
    class sequence_token
    {
    private:
        std::variant<tokens_t...> _selected_token{}; // we can not store the tokens here anymore!

    public:

        // need to set the token typerased?
        template <typename selected_token_t>
            requires (std::same_as<selected_token_t, tokens_t> || ...)
        sequence_token(selected_token_t && selected_token) :
            _selected_token{std::in_place_type<selected_token_t>, (selected_token_t &&)selected_token}
        {
        }
        sequence_token(sequence_token &&) = default;

        template <typename selected_token_t>
        sequence_token(selected_token_t const &) = delete;

    private:

        template <typename record_t>
        friend void tag_invoke(tag_t<libio::detokenize_to>, sequence_token &me, record_t &record)
        {
            std::visit([&] (auto & token) { libio::detokenize_to(token, record); }, me._selected_token);
        }
    };
} // namespace libio
