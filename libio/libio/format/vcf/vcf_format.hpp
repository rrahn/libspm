// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of the vcf format.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iosfwd>

#include <libio/format/vcf/vcf_token_header.hpp>
#include <libio/format/vcf/vcf_token_sample.hpp>
#include <libio/format/vcf/vcf_token_meta.hpp>
#include <libio/format/vcf/vcf_token.hpp>
#include <libio/format/vcf/vcf_header.hpp>
#include <libio/format/format_concept.hpp>
#include <libio/format/format_extension.hpp>
#include <libio/utility/tag_invoke.hpp>

namespace libio
{
    class vcf_format : public format_extension
    {
    private:

        vcf_header _header{};
    public:

        vcf_format() : format_extension{".vcf"}
        {
        }

        vcf_header const & header() const noexcept
        {
            return _header;
        }

    private:
        template <typename stream_t>
        constexpr friend auto tag_invoke(tag_t<libio::get_meta_token>, vcf_format const &, stream_t & istream) noexcept
            -> vcf_token_header<vcf_token_meta<stream_t>, vcf_token_sample<stream_t>>
        {
            return vcf_token_header{vcf_token_meta{istream}, vcf_token_sample{istream}};
        }

        template <typename stream_t>
        constexpr friend auto tag_invoke(tag_t<libio::format_token>, vcf_format const &, stream_t & istream)
        {
            return vcf_token{istream};
        }

        template <auto field_tag, typename ibuffer_t>
        constexpr friend auto tag_invoke(tag_t<libio::set_field>,
                                         vcf_format & me,
                                         field_code_type<field_tag> const & fc,
                                         ibuffer_t && ibuffer)
        {
            return libio::set_field(me._header, fc, (ibuffer_t &&)ibuffer);
        }
    };
} // namespace libio
