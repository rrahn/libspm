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

class vcf_parser_test : public ::testing::Test
{
public:
    std::filesystem::path reference_file{DATADIR"sim_ref_10Kb.fasta.gz"};

    auto load_sequences(std::filesystem::path reference_file, std::filesystem::path sample_sequence_file)
    {
        return std::pair{jstmap::load_sequences(reference_file).front(), jstmap::load_sequences(sample_sequence_file)};
    }

    template <typename jst_t, typename haplotypes_t>
    void test_jst_is_valid(jst_t const & jst, haplotypes_t const & haplotypes)
    {
        EXPECT_EQ(jst.size(), haplotypes.size());

        for (size_t haplotype_idx = 0; haplotype_idx < jst.size(); ++haplotype_idx)
            EXPECT_RANGE_EQ(jst.sequence_at(haplotype_idx), haplotypes[haplotype_idx]);
    }
};

TEST_F(vcf_parser_test, snps_only)
{
    std::filesystem::path vcf_file{DATADIR"sim_ref_10Kb_SNPs.vcf"};
    std::filesystem::path haplotype_file{DATADIR"sim_ref_10Kb_SNPs_haplotypes.fasta.gz"};

    jstmap::jst_t jst = std::move(jstmap::construct_jst_from_vcf(reference_file, vcf_file).front());

    auto [reference, haplotypes] = load_sequences(reference_file, haplotype_file);
    EXPECT_RANGE_EQ(jst.reference(), reference);
    test_jst_is_valid(jst, haplotypes);
}

TEST_F(vcf_parser_test, snps_and_indels)
{
    std::filesystem::path vcf_file{DATADIR"sim_ref_10Kb_SNP_INDELs.vcf"};
    std::filesystem::path haplotype_file{DATADIR"sim_ref_10Kb_SNP_INDELs_haplotypes.fasta.gz"};

    jstmap::jst_t jst = std::move(jstmap::construct_jst_from_vcf(reference_file, vcf_file).front());

    auto [reference, haplotypes] = load_sequences(reference_file, haplotype_file);
    EXPECT_RANGE_EQ(jst.reference(), reference);
    test_jst_is_valid(jst, haplotypes);
}

TEST_F(vcf_parser_test, sample_given_but_no_vcf_record)
{
    testing::internal::CaptureStderr();
    std::filesystem::path vcf_file{DATADIR"sim_ref_10Kb_no_variants.vcf"};

    EXPECT_TRUE(jstmap::construct_jst_from_vcf(reference_file, vcf_file).empty());

    auto expected_err_output = testing::internal::GetCapturedStderr();
    EXPECT_TRUE(expected_err_output.starts_with("[WARNING]"));
}

TEST_F(vcf_parser_test, unkown_reference_id)
{
    std::filesystem::path reference_file{DATADIR"in.fasta"};
    std::filesystem::path vcf_file{DATADIR"sim_ref_10Kb_SNP_INDELs.vcf"};

    EXPECT_THROW(jstmap::construct_jst_from_vcf(reference_file, vcf_file), std::runtime_error);
}

TEST_F(vcf_parser_test, empty_reference_file)
{
    std::filesystem::path reference_file{DATADIR"empty.fasta"};
    std::filesystem::path vcf_file{DATADIR"sim_ref_10Kb_SNP_INDELs.vcf"};

    EXPECT_THROW(jstmap::construct_jst_from_vcf(reference_file, vcf_file), std::runtime_error);
}
