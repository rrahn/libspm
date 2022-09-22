// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides vcf header.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <string>

#include <libio/file/field_code.hpp>
#include <libio/format/vcf/vcf_field_code.hpp>
namespace libio
{
    class vcf_header
    {
    private:
        std::string _version{};
        std::string _meta_infos{};
        std::string _sample_names{};
        size_t _sample_count{};

    public:
        constexpr vcf_header() = default;

        std::string const & version() const noexcept
        {
            return _version;
        }

        std::string const & infos() const noexcept
        {
            return _meta_infos;
        }

        std::string const & sample_names() const noexcept
        {
            std::cout << "Sample count: " << _sample_count << "\n";
            return _sample_names;
        }

    private:

        template <vcf_meta_field field_tag, typename chunked_buffer_t>
        friend auto tag_invoke(tag_t<libio::set_field>,
                               vcf_header &me,
                               field_code_type<field_tag>,
                               chunked_buffer_t &&buffer)
        {
            auto read_buffer = [&] (auto & target)
            {
                libio::read_token(target, buffer);
            };

            switch(field_tag) {
                case vcf_meta_field::version: read_buffer(me._version); break;
                case vcf_meta_field::meta: read_buffer(me._meta_infos); break;
                case vcf_meta_field::sample_names: read_buffer(me._sample_names); ++me._sample_count; break;
                default: break;//unkown field error!
            };
        }
    };
}  // namespace libio
