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
#include <libio/file/consume_tokenizer.hpp>
#include <libio/file/line_tokenizer.hpp>
#include <libio/file/segment_tokenizer.hpp>
#include <libio/format/fasta/fasta_field_code.hpp>
#include <libio/stream_token.hpp>
#include <libio/record/record_concept.hpp>

namespace libio
{
    template <typename record_t, typename token_t>
    concept fasta_record_c = requires (record_t & record, token_t & token)
    {
        libio::set_field(record, field_code<fasta_field::id>, token.get_area());
        libio::set_field(record, field_code<fasta_field::seq>, token.get_area());
    };

    template <typename stream_t>
    class fasta_token : protected stream_token<stream_t>
    {
    private:
        using base_t = stream_token<stream_t>;
        static constexpr auto token_signal_fn = seqan3::is_char<'>'>;
        static constexpr auto fragment_fn = seqan3::is_char<'\n'> || seqan3::is_char<'\r'>;

    public:
        // here we can describe whether we want to skip
        fasta_token(stream_t &stream) : base_t{stream, token_signal_fn}
        {
        }
        fasta_token(fasta_token const &) = delete;
        fasta_token(fasta_token &&) = default;

    private:
        // what does this mean? how can materialise this now?
        template <typename record_t>
            requires fasta_record_c<record_t, fasta_token>
        friend void tag_invoke(tag_t<libio::detokenize_to>, fasta_token &me, record_t &record)
        {
            libio::set_field(record, field_code<fasta_field::id>, consume_tokenizer{line_tokenizer{me.get_area()}});
            libio::set_field(record,
                             field_code<fasta_field::seq>,
                             consume_tokenizer{segment_tokenizer{me.get_area(), seqan3::is_alpha}});
        }
    };

    template <typename stream_buffer_t>
    fasta_token(stream_buffer_t *) -> fasta_token<stream_buffer_t>;
} // namespace libio
