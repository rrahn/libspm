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
#include <libio/file/until_tokenizer.hpp>
#include <libio/format/vcf/vcf_field_code.hpp>
#include <libio/stream_token.hpp>
#include <libio/record/record_concept.hpp>

namespace libio
{
    // maybe line_token is different -> because line-token has signal at end.
    // using vcf_token_tag_t  = decltype(until_token{seqan3::is_char<'\n'> || seqan3::is_char<'\r'>});

    template <typename stream_t>
    class vcf_token : protected stream_token<stream_t, line_token>
    {
    private:
        using base_t = stream_token<stream_t, line_token>;

        static constexpr auto delimiter_fn = seqan3::is_char<'\t'>;
    public:
        // here we can describe whether we want to skip
        vcf_token(stream_t &stream) : base_t{stream, line_token{}}
        {
        }
        vcf_token(vcf_token const &) = delete;
        vcf_token(vcf_token &&) = default;

    private:
        // what does this mean? how can materialise this now?
        template <typename record_t>
        friend void tag_invoke(tag_t<libio::detokenize_to>, vcf_token &me, record_t &record)
        {
            auto next_field = [&] ()
            {
                return consume_tokenizer{until_tokenizer{me.get_area(), delimiter_fn}};
            };

            [&] <size_t ...idx> (std::index_sequence<idx...> const &)
            {
                (libio::set_field(record, field_code<static_cast<vcf_field>(idx + 1)>, next_field()),...);
            } (std::make_index_sequence<static_cast<int32_t>(vcf_field::genotype_format)>());

            // we need to read the last field!
            while (me.get_area().begin() != me.get_area().end())
                libio::set_field(record,
                                 field_code<vcf_field::genotypes>,
                                 consume_tokenizer{until_tokenizer{me.get_area(), delimiter_fn}});
        }
    };

    template <typename stream_t>
    vcf_token(stream_t &) -> vcf_token<stream_t>;
} // namespace libio
