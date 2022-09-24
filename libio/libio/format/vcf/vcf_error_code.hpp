// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides error code adaption for the vcf format.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <system_error>

#include <libio/format/vcf/vcf_field_code.hpp>

namespace libio
{
    // ----------------------------------------------------------------------------
    // VCF meta field error defintion
    // ----------------------------------------------------------------------------

    namespace _vcf_meta_field_error
    {
        class category final : public std::error_category
        {
        public:
            category() = default;

            virtual const char *name() const noexcept override
            {
                return "vcf_meta_field";
            }

            virtual std::string message(int ev) const override
            {
                using namespace std::literals;
                switch(static_cast<vcf_meta_field>(ev))
                {
                    case vcf_meta_field::version: return "no valid version field."s;
                    case vcf_meta_field::meta: return "no valid meta field."s;
                    case vcf_meta_field::sample_names: return "no valid sample names field."s;
                    default: return "[unknown error]"s;
                }
            }
        };
    } // namespace _vcf_meta_field_error

    inline const _vcf_meta_field_error::category vcf_meta_field_error_category{};

    // Customise make_error CPO!
    std::error_code make_error_code(vcf_meta_field error)
    {
        return {static_cast<int>(error), vcf_meta_field_error_category};
    }

    // ----------------------------------------------------------------------------
    // VCF field error defintion
    // ----------------------------------------------------------------------------

    namespace _vcf_field_error
    {
        class category final : public std::error_category
        {
        public:
            category() = default;

            virtual const char *name() const noexcept override
            {
                return "vcf_field";
            }

            virtual std::string message(int ev) const override
            {
                using namespace std::literals;
                switch(static_cast<vcf_field>(ev))
                {
                    case vcf_field::chrom: return "no valid chromosome field."s;
                    case vcf_field::pos: return "no valid pos field."s;
                    case vcf_field::id: return "no valid id field."s;
                    case vcf_field::ref: return "no valid ref field."s;
                    case vcf_field::alt: return "no valid alt field."s;
                    case vcf_field::qual: return "no valid qual field."s;
                    case vcf_field::filter: return "no valid filter field."s;
                    case vcf_field::info: return "no valid info field."s;
                    case vcf_field::genotype_format: return "no valid format field."s;
                    case vcf_field::genotypes: return "no valid genotypes field."s;
                    default: return "[unknown error]"s;
                }
            }
        };
    } // namespace _vcf_field_error

    inline const _vcf_field_error::category vcf_field_error_category{};

    // Customise make_error CPO!
    std::error_code make_error_code(vcf_field error)
    {
        return {static_cast<int>(error), vcf_field_error_category};
    }
}  // namespace libio

namespace std
{
    template <>
    struct is_error_code_enum<libio::vcf_meta_field> : true_type {};

    template <>
    struct is_error_code_enum<libio::vcf_field> : true_type {};
} // namespace std
