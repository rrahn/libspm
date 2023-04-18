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

#include <libjst/coverage/range_domain.hpp>

#include <jstmap/global/jstmap_types.hpp>

namespace jstmap
{
    class stripped_vcf_record
    {
    private:
        using alternative_t = std::string;
        using genotypes_t = std::vector<coverage_t>;
        using position_t = uint32_t;
        using coverage_value_t = typename coverage_t::value_type;
        using domain_t = libjst::range_domain<coverage_value_t>;

        seqan::VcfIOContext<> const *_io_context{};

        // relevant fields to set.
        std::string _ref{};
        std::vector<alternative_t> _alt{};
        genotypes_t _genotypes{};
        domain_t _domain{};
        position_t _pos{};
        int16_t _chrom{};
        uint8_t _alternative_count{};

    public:

        using intermediate_variant_t = std::tuple<int32_t, std::string, int32_t, coverage_t>;

        stripped_vcf_record() = default;
        stripped_vcf_record(seqan::VcfRecord &,
                            seqan::VcfIOContext<> const &,
                            domain_t);

        std::string_view contig_name() const noexcept;
        genotypes_t const & field_genotype() const noexcept;
        void alternatives(rcs_store_t &) const;

    private:
        void set_field_chrom(seqan::VcfRecord &) noexcept;

        void set_field_pos(seqan::VcfRecord &) noexcept;

        void set_field_ref(seqan::VcfRecord &);

        void set_field_alt(seqan::VcfRecord &);

        void set_field_genotype(seqan::VcfRecord &);
    };
}  // namespace jstmap
