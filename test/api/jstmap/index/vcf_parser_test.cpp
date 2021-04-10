// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

// Including required debug stream types for range_eq test.
#include <seqan3/core/detail/debug_stream_alphabet.hpp>
#include <seqan3/core/detail/debug_stream_range.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include <jstmap/index/load_sequence.hpp>
#include <jstmap/index/vcf_parser.hpp>

TEST(vcf_parser_test, snps_only)
{
    std::filesystem::path reference_file{DATADIR"sim_ref_10Kb.fasta.gz"};
    std::filesystem::path vcf_file{DATADIR"sim_ref_10Kb_SNPs.vcf"};
    std::filesystem::path sample_sequence_file{DATADIR"sim_ref_10Kb_SNPs_haplotypes.fasta.gz"};

    std::vector haplotypes = jstmap::load_sequences(sample_sequence_file);
    std::vector reference = jstmap::load_sequences(reference_file);
    jstmap::jst_t jst = jstmap::construct_jst_from_vcf(reference_file, vcf_file);

    EXPECT_EQ(haplotypes.size(), jst.size());
    EXPECT_RANGE_EQ(jst.reference(), reference[0]);

    for (size_t haplotype_idx = 0; haplotype_idx < jst.size(); ++haplotype_idx)
        EXPECT_RANGE_EQ(jst.sequence_at(haplotype_idx), haplotypes[haplotype_idx]);
}
