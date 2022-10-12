// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the function to parse a vcf file and construct a JST from it.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <utility>
#include <vector>
#include <tuple>

#include <seqan/vcf_io.h>

#include <libjst/utility/bit_vector.hpp>

#include <jstmap/global/jstmap_jst_types.hpp>

namespace jstmap
{
    class stripped_vcf_record
    {
    private:
        using alternative_t = std::string;
        using coverage_t = libjst::bit_vector<>;
        using genotypes_t = std::vector<coverage_t>;
        using position_t = libjst::variant_position_t<snp_t>;

        seqan::VcfIOContext<> const *_io_context{};

        // relevant fields to set.
        std::string _ref{};
        std::vector<alternative_t> _alt{};
        genotypes_t _genotypes{};
        position_t _pos{};
        int16_t _chrom{};
        uint8_t _alternative_count{};

    public:

        using intermediate_variant_t = std::tuple<int32_t, std::string, int32_t, coverage_t>;

        stripped_vcf_record() = default;
        stripped_vcf_record(seqan::VcfRecord &, seqan::VcfIOContext<> const &);

        std::string_view contig_name() const noexcept;
        genotypes_t const & field_genotype() const noexcept;
        void alternatives(variant_store_t &) const;

    private:
        void set_field_chrom(seqan::VcfRecord &) noexcept;

        void set_field_pos(seqan::VcfRecord &) noexcept;

        void set_field_ref(seqan::VcfRecord &);

        void set_field_alt(seqan::VcfRecord &);

        void set_field_genotype(seqan::VcfRecord &);
    };
}  // namespace jstmap
