// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <filesystem>

#include <jstmap/view/load_jst.hpp>
#include <jstmap/view/view_format_fasta.hpp>

TEST(jstmap_view_test, view_format_fasta)
{
    std::filesystem::path jst_file{DATADIR"sim_refx5.jst"};
    jstmap::jst_t jst = jstmap::load_jst(jst_file);

    testing::internal::CaptureStdout();

    jstmap::view_as_format(jst, 0);

    std::string captured_output = testing::internal::GetCapturedStdout();

    std::string expected_output{R"fasta(> ID_0
TATGCACCAGAGTATGGAAGCATAAGCTCTGCATGCAAAGGTACATCAGATCCTGCGGTTGGGTGCCAACCCAAGTGTGT
TCACGGGCGC
)fasta"};

    EXPECT_EQ(captured_output, expected_output);
}

TEST(jstmap_view_test, view_format_fasta_unknown_haplpotype_index)
{
    std::filesystem::path jst_file{DATADIR"sim_refx5.jst"};
    jstmap::jst_t jst = jstmap::load_jst(jst_file);

    EXPECT_THROW(jstmap::view_as_format(jst, 6), std::out_of_range);
}
