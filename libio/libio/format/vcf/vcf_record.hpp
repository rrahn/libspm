// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides vcf record.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <charconv>
#include <string>

#include <libio/file/field_code.hpp>
#include <libio/format/vcf/vcf_field_code.hpp>
namespace libio
{
    class vcf_record
    {
    private:
        // std::string _chrom{};
        // std::string _pos{};
        // std::string _id{};
        // std::string _ref{};
        // std::string _alt{};
        // std::string _qual{};
        // std::string _filter{};
        // std::string _info{};
        // std::string _genotype_format{};
        // std::string _genotypes{};

        std::string _values{};
        std::vector<size_t> _offsets{};

    public:
        vcf_record() noexcept
        {
            _offsets.resize(11, 0);
        }

        int32_t chrom() const noexcept
        {
            int32_t _number{};
            std::string_view _chrom = get_field(field_code<vcf_field::chrom>);
            std::from_chars(_chrom.data(), _chrom.data() + _chrom.size(), _number);
            return _number;
        }

        std::string_view pos() const noexcept
        {
            return get_field(field_code<vcf_field::pos>);
        }

        std::string_view id() const noexcept
        {
            return get_field(field_code<vcf_field::id>);
        }

        std::string_view ref() const noexcept
        {
            return get_field(field_code<vcf_field::ref>);
        }

        std::string_view alt() const noexcept
        {
            return get_field(field_code<vcf_field::alt>);
        }

        std::string_view qual() const noexcept
        {
            return get_field(field_code<vcf_field::qual>);
        }

        std::string_view filter() const noexcept
        {
            return get_field(field_code<vcf_field::filter>);
        }

        std::string_view info() const noexcept
        {
            return get_field(field_code<vcf_field::info>);
        }

        std::string_view genotype_format() const noexcept
        {
            return get_field(field_code<vcf_field::genotype_format>);
        }

        std::string_view genotypes() const noexcept
        {
            return get_field(field_code<vcf_field::genotypes>);
        }

        void clear() noexcept
        {
            _values.clear();
            std::ranges::for_each(_offsets, [] (size_t & o) { o = 0; });
        }

        size_t bytes() const noexcept
        {
            return _values.size();
        }
    private:

        template <vcf_field field_tag>
        std::string_view get_field(field_code_type<field_tag> const &fc) const noexcept
        {
            return {_values.data() + _offsets[fc.index()-1], _values.data() + _offsets[fc.index()]};
        }

        template <vcf_field field_tag, typename chunked_buffer_t>
        friend auto tag_invoke(tag_t<libio::set_field>,
                               vcf_record &me,
                               field_code_type<field_tag> const &fc,
                               chunked_buffer_t &&buffer)
        {
            // auto read_buffer = [&] (auto & target)
            // {
                libio::read_token(me._values, buffer);
                me._offsets[fc.index()] = me._values.size();
            // };

            // switch(field_tag) {
            //     case vcf_field::chrom: read_buffer(me._chrom); break;
            //     case vcf_field::pos: read_buffer(me._pos); break;
            //     case vcf_field::id: read_buffer(me._id); break;
            //     case vcf_field::ref: read_buffer(me._ref); break;
            //     case vcf_field::alt: read_buffer(me._alt); break;
            //     case vcf_field::qual: read_buffer(me._qual); break;
            //     case vcf_field::filter: read_buffer(me._filter); break;
            //     case vcf_field::info: read_buffer(me._info); break;
            //     case vcf_field::genotype_format: read_buffer(me._genotype_format); break;
            //     case vcf_field::genotypes: read_buffer(me._genotypes); break;
            //     default: break;//unkown field error!
            // };
        }
    };
}  // namespace libio
