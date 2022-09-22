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

#include <seqan3/utility/char_operations/predicate.hpp>

#include <libio/file/field_code.hpp>
#include <libio/file/consume_tokenizer.hpp>
#include <libio/file/line_tokenizer.hpp>
#include <libio/file/segment_tokenizer.hpp>
#include <libio/file/until_tokenizer.hpp>
#include <libio/format/fastq/fastq_field_code.hpp>
#include <libio/record/record_concept.hpp>
#include <libio/stream_token.hpp>

namespace libio
{
    template <typename record_t, typename token_t>
    concept fastq_record_c = requires (record_t & record, token_t & token)
    {
        libio::set_field(record, field_code<fastq_field::id>, token.get_area());
        libio::set_field(record, field_code<fastq_field::seq>, token.get_area());
        libio::set_field(record, field_code<fastq_field::qual>, token.get_area());
    };


    using fastq_token_tag_t = decltype(until_token{seqan3::is_char<'@'>});

    template <typename stream_t>
    class fastq_token : protected stream_token<stream_t, fastq_token_tag_t>
    {
    private:

        using base_t = stream_token<stream_t, fastq_token_tag_t>;
        static constexpr auto is_qual_signal = seqan3::is_char<'+'>;
        static constexpr auto is_newline = seqan3::is_char<'+'>;
        static constexpr auto is_phred = seqan3::is_in_interval<33, 126>;

    public:
        // here we can describe whether we want to skip
        explicit fastq_token(stream_t &stream) : base_t{stream, fastq_token_tag_t{}}
        {
        }
        fastq_token(fastq_token const &) = delete;
        fastq_token(fastq_token &&) = default;

    private:
        // what does this mean? how can materialise this now?
        template <typename record_t>
            requires fastq_record_c<record_t, fastq_token>
        friend auto tag_invoke(tag_t<libio::detokenize_to>, fastq_token &me, record_t &record)
        {
            libio::set_field(record, field_code<fastq_field::id>, consume_tokenizer{line_tokenizer{me.get_area()}});
            libio::set_field(record, field_code<fastq_field::seq>, consume_tokenizer{line_tokenizer{me.get_area()}});
            // Make an assert tokenizer
            std::ranges::for_each(line_tokenizer{me.get_area()},
                                  [] ([[maybe_unused]] auto && seq) { assert(seq.empty() || is_qual_signal(seq[0])); });
            libio::set_field(record, field_code<fastq_field::qual>, consume_tokenizer{segment_tokenizer{me.get_area(), is_phred}});
        }
    };

    template <typename stream_t>
    fastq_token(stream_t &) -> fastq_token<stream_t>;
} // namespace libio
