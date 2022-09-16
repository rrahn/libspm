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
#include <libio/file/line_tokenizer.hpp>
#include <libio/file/segment_tokenizer.hpp>
#include <libio/file/until_tokenizer.hpp>
#include <libio/format/fastq/fastq_field_code.hpp>
#include <libio/record/record_concept.hpp>
#include <libio/stream_token.hpp>

namespace libio
{
    template <typename streambuffer_t>
    class fastq_token : protected stream_token<streambuffer_t>
    {
    private:

        using base_t = stream_token<streambuffer_t>;
        static constexpr auto is_token_signal = seqan3::is_char<'@'>;
        static constexpr auto is_qual_signal = seqan3::is_char<'+'>;
        static constexpr auto is_newline = seqan3::is_char<'+'>;
        static constexpr auto is_phred = seqan3::is_in_interval<33, 126>;

    public:
        // here we can describe whether we want to skip
        fastq_token(streambuffer_t *stream_buffer) : base_t{stream_buffer, is_token_signal}
        {
        }

    private:
        // what does this mean? how can materialise this now?
        template <typename record_t>
        friend auto tag_invoke(tag_t<libio::detokenize_to>, fastq_token &me, record_t &record)
        {
            libio::set_field(record, field_code<fastq_field::id>, line_tokenizer{me.get_area()});
            libio::set_field(record, field_code<fastq_field::seq>, line_tokenizer{me.get_area()});
            // Make an assert tokenizer
            std::ranges::for_each(line_tokenizer{me.get_area()},
                                  [] ([[maybe_unused]] auto && seq) { assert(seq.empty() || is_qual_signal(seq[0])); });
            libio::set_field(record, field_code<fastq_field::qual>, segment_tokenizer{me.get_area(), is_phred});
        }
    };

    template <typename stream_buffer_t>
    fastq_token(stream_buffer_t *) -> fastq_token<stream_buffer_t>;
} // namespace libio
