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

    void stripped_vcf_record::alternatives(rcs_store_t & store) const
    {
        if (_alternative_count != _genotypes.size())
            throw std::logic_error{"Invalid number of coverages and alternative count."};

        for (size_t i = 0; i < _alternative_count; ++i)
        {
            if (_alt[i][0] == '<') continue; // skip these alternatives for now

            reference_t alt_sequence{_alt[i].begin(), _alt[i].end()};
            variant_t variant{};

            if (_ref.size() == _alt[i].size() && _ref.size() == 1) { // SNP
                ++_stat->snv_count;
                variant = variant_t{libjst::breakpoint{_pos, 1}, std::move(alt_sequence), std::move(_genotypes[i])};
            } else { // generic alternative: SNP, InDel, etc.
                ++_stat->indel_count;
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
            if (!store.variants().has_conflicts(variant)) {
                store.add(std::move(variant));
            } else {
                std::cout << "CONFLICT\n";
            }
        }
    }

    stripped_vcf_record::genotypes_t const & stripped_vcf_record::field_genotype() const noexcept
    {
        return _genotypes;
    }

    void stripped_vcf_record::set_field_chrom(std::string_view field)
    {
        _chrom_name.assign(field.begin(), field.end());
    }

    void stripped_vcf_record::set_field_pos(std::string_view field)
    {
        if (auto res = std::from_chars(field.data(), field.data() + field.size(), _pos); res.ptr == field.data())
            throw std::make_error_code(res.ec);

        --_pos; // correct 1-based offset.
    }

    void stripped_vcf_record::set_field_ref(std::string_view field)
    {
        _ref.assign(field.begin(), field.end());
    }

    void stripped_vcf_record::set_field_alt(std::string_view field)
    {
        _alternative_count = 1;
        auto first = field.begin();
        for (auto it = field.begin(); it != field.end(); ++it) { // go over alternative
            if (*it == ',') { // found a new alternative
                _alt.emplace_back(first, it); // emplace the alternative
                ++_alternative_count;
                ++it; // jump over comma
                first = it; // set next begin.
            }
        }
        // emplace last element
        assert(first != field.end());
        _alt.emplace_back(first, field.end()); // emplace the alternative
    }

    void stripped_vcf_record::set_field_genotype(std::string_view genotypes)
    {
        _genotypes.resize(_alternative_count, coverage_t{_domain});

        auto record_coverage = [&] (auto alt_first, auto alt_last, size_t const haplotype_idx) {
            assert(haplotype_idx < _haplotype_count);
            int32_t alt_index{};
            if (std::from_chars(std::to_address(alt_first), std::to_address(alt_last), alt_index).ec != std::errc{})
                throw std::runtime_error{"Extracting haplotype failed!"};

            if (alt_index > 0) {
                coverage_t & current_coverage = _genotypes[--alt_index];
                current_coverage.insert(current_coverage.end(), haplotype_idx);
            }
        };

        size_t haplotype_idx{};
        for (std::ptrdiff_t sample_idx = 0; sample_idx < _sample_count; ++sample_idx) {
            // TODO: determine ploidy
            std::string_view genotype = read_field(genotypes);
            auto genotype_it = genotype.begin();
            // TODO: fix if more than 9 alternatives per position.
            record_coverage(genotype_it, genotype_it + 1, haplotype_idx++);
            genotype_it += 2;
            record_coverage(genotype_it, genotype_it + 1, haplotype_idx++);
        }
    }

    std::string_view stripped_vcf_record::read_field(std::string_view & buffer) noexcept {
        auto delimiter_ptr = std::memchr(std::to_address(buffer.begin()), '\t', buffer.size());
        if (delimiter_ptr == nullptr) {
            std::string_view field{};
            std::swap(buffer, field);
            return field;
        }

        char const * field_end = reinterpret_cast<char const *>(delimiter_ptr);
        std::ptrdiff_t field_span = field_end - buffer.data();
        std::string_view field{buffer.begin(), buffer.begin() + field_span};
        buffer = buffer.substr(field_span + 1);
        return field;
    }

}  // namespace jstmap
