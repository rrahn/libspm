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

#include <algorithm>
#include <charconv>
#include <iostream>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>

#include <jstmap/create/stripped_vcf_record.hpp>

namespace jstmap
{
    stripped_vcf_record::stripped_vcf_record(seqan::VcfRecord & record,
                                             seqan::VcfIOContext<> const & io_context) :
        _io_context{std::addressof(io_context)}
    {
        set_field_chrom(record);
        set_field_pos(record);
        set_field_ref(record);
        set_field_alt(record);
        set_field_genotype(record);
    }

    void stripped_vcf_record::alternatives(rcs_store_t & store) const
    {
        if (_alternative_count != _genotypes.size())
            throw std::logic_error{"Invalid number of coverages and alternative count."};

        for (size_t i = 0; i < _alternative_count; ++i)
        {
            if (_alt[i][0] == '<') continue; // skip these alternatives for now
            if (_ref.size() == _alt[i].size() && _ref.size() == 1) { // SNP
                store.add(_pos, variant_t{alphabet_t{_alt[i][0]}}, std::move(_genotypes[i]));
                // store.emplace(snp_t{_pos, }, std::move(_genotypes[i]));
            } else { // generic alternative: SNP, InDel, etc.
                continue; // skip for now!
                // auto [fst_ref, fst_alt] = std::ranges::mismatch(_ref, _alt[i]); // first non equal ranges in the beginning
                // auto ref_suffix = _ref | std::views::drop(std::ranges::distance(std::ranges::begin(_ref), fst_ref))
                //                        | std::views::reverse;
                // auto alt_suffix = _alt[i] | std::views::drop(std::ranges::distance(std::ranges::begin(_alt[i]), fst_alt))
                //                           | std::views::reverse;
                // auto [lst_ref_rev, lst_alt_rev] = std::ranges::mismatch(ref_suffix, alt_suffix);

                // auto allele_view = std::ranges::subrange{fst_alt, lst_alt_rev.base()} | std::views::transform(
                //     [] (char c) { return alphabet_t{c}; });
                // std::vector<alphabet_t> allele{std::ranges::begin(allele_view), std::ranges::end(allele_view)};
                // libjst::variant_deletion_t<indel_t> deletion = std::ranges::distance(fst_ref, lst_ref_rev.base());
                // store.emplace(indel_t{_pos, std::move(allele), deletion}, std::move(_genotypes[i]));
            }
        }
    }

    stripped_vcf_record::genotypes_t const & stripped_vcf_record::field_genotype() const noexcept
    {
        return _genotypes;
    }

    std::string_view stripped_vcf_record::contig_name() const noexcept
    {
        assert(_io_context != nullptr);
        assert(_chrom < seqan::length(seqan::contigNames(*_io_context)));
        return {seqan::toCString(seqan::contigNames(*_io_context)[_chrom])};
    }

    void stripped_vcf_record::set_field_chrom(seqan::VcfRecord & record) noexcept
    {
        _chrom = record.rID;
    }

    void stripped_vcf_record::set_field_pos(seqan::VcfRecord & record) noexcept
    {
        _pos = record.beginPos;
    }

    void stripped_vcf_record::set_field_ref(seqan::VcfRecord & record)
    {
        _ref.assign(seqan::begin(record.ref, seqan::Standard{}), seqan::end(record.ref, seqan::Standard{}));
    }

    void stripped_vcf_record::set_field_alt(seqan::VcfRecord & record)
    {
        // using std::swap;
        // swap(_alt, record.alt);

        _alternative_count = 1;
        auto first = seqan::begin(record.alt);
        for (auto it = seqan::begin(record.alt); it != seqan::end(record.alt); ++it) { // go over alternative
            if (*it == ',') { // found a new alternative
                _alt.emplace_back(first, it); // emplace the alternative
                ++_alternative_count;
                ++it; // jump over comma
                first = it; // set next begin.
            }
        }
        // emplace last element
        assert(first != seqan::end(record.alt));
        _alt.emplace_back(first, seqan::end(record.alt)); // emplace the alternative
    }

    void stripped_vcf_record::set_field_genotype(seqan::VcfRecord & record)
    {
        size_t total_haplotypes = seqan::length(seqan::sampleNames(*_io_context)) << 1;
        coverage_t coverage{};
        coverage.resize(total_haplotypes, false);
        _genotypes.resize(_alternative_count, std::move(coverage));

        auto set_coverage = [&] (char const * ptr, size_t const haplotype_count) {
            assert(haplotype_count < total_haplotypes);
            uint32_t alt_index{};
            if (std::from_chars(ptr, ptr + 1, alt_index).ec != std::errc{})
                throw std::runtime_error{"Extracting haplotype failed!"};

            if (alt_index > 0)
                _genotypes[--alt_index][haplotype_count] = true;
        };

        size_t haplotype_count{};
        for (auto it = seqan::begin(record.genotypeInfos); it != seqan::end(record.genotypeInfos); ++it)
        {
            set_coverage(std::addressof((*it)[0]), haplotype_count++);
            set_coverage(std::addressof((*it)[2]), haplotype_count++);
        }
    }
}  // namespace jstmap
