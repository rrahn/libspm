// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implements the token for the fasta format.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/utility/char_operations/predicate.hpp>

#include <libio/file/field_code.hpp>
#include <libio/file/line_tokenizer.hpp>
#include <libio/file/segment_tokenizer.hpp>
#include <libio/format/fasta/fasta_field_code.hpp>
#include <libio/stream_token.hpp>
#include <libio/record/record_concept.hpp>

namespace libio
{
    template <typename streambuffer_t>
    class fasta_token : protected stream_token<streambuffer_t>
    {
    private:
        using base_t = stream_token<streambuffer_t>;
        static constexpr auto token_signal_fn = seqan3::is_char<'>'>;
        static constexpr auto fragment_fn = seqan3::is_char<'\n'> || seqan3::is_char<'\r'>;

    public:
        // here we can describe whether we want to skip
        fasta_token(streambuffer_t *stream_buffer) : base_t{stream_buffer, token_signal_fn}
        {
        }

    private:
        // what does this mean? how can materialise this now?
        template <typename record_t>
        friend auto tag_invoke(tag_t<libio::detokenize_to>, fasta_token &me, record_t &record)
        {
            libio::set_field(record, field_code<fasta_field::id>, line_tokenizer{me.get_area()});
            libio::set_field(record, field_code<fasta_field::seq>, segment_tokenizer{me.get_area(), seqan3::is_alpha});
        }
    };

    template <typename stream_buffer_t>
    fasta_token(stream_buffer_t *) -> fasta_token<stream_buffer_t>;
} // namespace libio
