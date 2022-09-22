// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of the vcf header structure.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/utility/char_operations/predicate.hpp>

#include <libio/file/consume_tokenizer.hpp>
#include <libio/file/line_tokenizer.hpp>
#include <libio/file/until_tokenizer.hpp>
#include <libio/file/tokenization.hpp>
#include <libio/format/vcf/vcf_field_code.hpp>
#include <libio/record/record_concept.hpp>
#include <libio/stream_token.hpp>
#include <libio/utility/tag_invoke.hpp>

namespace libio
{
    using vcf_token_sample_tag_t = decltype(until_token{seqan3::is_char<'\n'> || seqan3::is_char<'\r'>});

    template <typename stream_t>
    class vcf_token_sample : protected stream_token<stream_t, line_token>
    {
    private:
        using base_t = stream_token<stream_t, line_token>;
        // requires now some pattern matching
        static constexpr auto delimiter_fn = seqan3::is_char<'\t'>;
    public:
        // here we can describe whether we want to skip
        vcf_token_sample(stream_t &stream) : base_t{stream, line_token{}}
        {
        }
        vcf_token_sample(vcf_token_sample const &) = delete;
        vcf_token_sample(vcf_token_sample &&) = default;

        using base_t::get_area;
    private:
        // what does this mean? how can materialise this now?
        template <typename header_t>
        friend void tag_invoke(tag_t<libio::detokenize_to>, vcf_token_sample &me, header_t &header)
        {
            // we need to think about this.
            // split and count number of splits?
            for (int32_t i = 0; i < field_code_type<vcf_field::genotype_format>::index(); ++i)
                consume_tokenizer{until_tokenizer{me.get_area(), delimiter_fn}};
            // read sample names and call it until no more samples are available.
            // the header can then move them all into some other region.
            while (me.get_area().begin() != me.get_area().end())
                libio::set_field(header,
                                 field_code<vcf_meta_field::sample_names>,
                                 consume_tokenizer{until_tokenizer{me.get_area(), delimiter_fn}});
        }
    };

    template <typename stream_t>
    vcf_token_sample(stream_t &) -> vcf_token_sample<stream_t>;
} // namespace libio
