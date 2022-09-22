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

#include <libio/file/consume_tokenizer.hpp>
#include <libio/file/pivot_tokenizer.hpp>
#include <libio/file/line_tokenizer.hpp>
#include <libio/file/until_tokenizer.hpp>
#include <libio/file/tokenization.hpp>
#include <libio/format/vcf/vcf_field_code.hpp>
#include <libio/stream_token.hpp>
#include <libio/utility/tag_invoke.hpp>

namespace libio
{

    using vcf_meta_token_tag_t = decltype(pivot_token{pivot_matcher{"#CHROM"}});
    // we need another abstraction to set the token delimiter.
    template <typename stream_t>
    class vcf_token_meta : protected stream_token<stream_t, vcf_meta_token_tag_t>
    {
    private:
        using base_t = stream_token<stream_t, vcf_meta_token_tag_t>;

        static constexpr auto delimiter_fn = seqan3::is_char<'#'>;  // pattern: ##
    public:
        // here we can describe whether we want to skip
        vcf_token_meta(stream_t &stream) : base_t{stream, pivot_token{pivot_matcher{"#CHROM"}}}
        {
        }
        vcf_token_meta(vcf_token_meta const &) = delete;
        vcf_token_meta(vcf_token_meta &&) = default;

        using base_t::get_area;

    private:
        // what does this mean? how can materialise this now?
        template <typename header_t>
        friend void tag_invoke(tag_t<libio::detokenize_to>, vcf_token_meta &me, header_t &header)
        {
            // what if the next values not available?
            auto _set_field = [&] <vcf_meta_field field_tag> (field_code_type<field_tag> fc)
            {
                auto consume_field = consume_tokenizer{until_tokenizer{me.get_area(), delimiter_fn}};
                libio::set_field(header, fc, line_tokenizer{consume_field});
            };

            // skip beginning until found start_tag -> skip until find not start_tag -> read until end_tag
            _set_field(field_code<vcf_meta_field::version>);

            // get_area is cached so we consume initial area only once!
            while (me.get_area().begin() != me.get_area().end())
                _set_field(field_code<vcf_meta_field::meta>);
                // libio::set_field(header, field_code<vcf_meta_field::meta>, mata_field());
        }
    };

    template <typename stream_t>
    vcf_token_meta(stream_t &) -> vcf_token_meta<stream_t>;
} // namespace libio
