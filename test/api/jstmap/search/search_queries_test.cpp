// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <jstmap/global/load_jst.hpp>
#include <jstmap/search/load_queries.hpp>
#include <jstmap/search/search_queries.hpp>

// TEST(jstmap_index, search_jst)
// {
//     std::filesystem::path jst_file{DATADIR"sim_refx5.jst"};
//     jstmap::jst_t jst = jstmap::load_jst(jst_file);

//     std::filesystem::path queries_file{DATADIR"sim_reads_ref1x10.fa"};
//     std::vector reads = jstmap::load_queries(queries_file);

//     std::vector results = jstmap::search_queries(std::move(jst), std::move(reads));

//     EXPECT_TRUE((std::ranges::find(results, libjst::context_position{0, 00})) != results.end());
//     EXPECT_TRUE((std::ranges::find(results, libjst::context_position{0, 16})) != results.end());
//     EXPECT_TRUE((std::ranges::find(results, libjst::context_position{0, 36})) != results.end());
//     EXPECT_TRUE((std::ranges::find(results, libjst::context_position{0, 01})) != results.end());
//     EXPECT_TRUE((std::ranges::find(results, libjst::context_position{0, 21})) != results.end());
//     EXPECT_TRUE((std::ranges::find(results, libjst::context_position{0, 41})) != results.end());
//     EXPECT_TRUE((std::ranges::find(results, libjst::context_position{0, 61})) != results.end());
//     EXPECT_TRUE((std::ranges::find(results, libjst::context_position{0, 70})) != results.end());
//     EXPECT_TRUE((std::ranges::find(results, libjst::context_position{0, 41})) != results.end());
//     EXPECT_TRUE((std::ranges::find(results, libjst::context_position{0, 50})) != results.end());

//     EXPECT_EQ(results.size(), 10u);
// }

TEST(jstmap_index, search_jst_2)
{
    // std::filesystem::path jst_file{DATADIR"sim_refx5_p0.jst"};

    // auto [jst, partitioned_jst_handle] = jstmap::load_jst(jst_file);
    // jstmap::partitioned_jst_t const & partitioned_jst = *partitioned_jst_handle;

    // std::filesystem::path queries_file{DATADIR"sim_reads_ref1x10.fa"};
    // std::vector reads = jstmap::load_queries(queries_file);

    // std::vector results = jstmap::search_queries(std::move(jst), std::move(reads));

    using seqan3::operator""_dna5;

    using jst_t = libjst::journaled_sequence_tree<jstmap::raw_sequence_t>;
    using partitioned_jst_t = libjst::journal_sequence_tree_partitioned<jst_t>;
    // using event_t = typename jst_t::event_type;
    // using snp_t = typename event_t::snp_type;
    // using substitution_t = typename event_t::substitution_type;
    // using coverage_t = typename event_t::coverage_type;

                //                0         1         2         3         4         5         6         7         8         9
                //                0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    jstmap::raw_sequence_t ref = "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt"_dna5;
    //                  aa
    //                    c
    jst_t jst{std::move(ref), 4};
    partitioned_jst_t pjst{std::addressof(jst)};

    jstmap::raw_sequence_t q1{"cgtacgtacgtacgtacgtacgta"_dna5};
    jstmap::raw_sequence_t q2{"acgtacgtacgtacgtacgtacgt"_dna5};
    jstmap::raw_sequence_t q3{"gtacgtacgtacgtacgtacgtac"_dna5};

    seqan::StringSet<std::views::all_t<jstmap::raw_sequence_t const &>> pattern_collection{};
    seqan::appendValue(pattern_collection, std::as_const(q1) | std::views::all);
    seqan::appendValue(pattern_collection, std::as_const(q2) | std::views::all);
    seqan::appendValue(pattern_collection, std::as_const(q3) | std::views::all);

    // jst.insert(event_t{2ull, substitution_t{"AA"_dna5}, coverage_t{1, 0, 1, 0}});
    // jst.insert(event_t{4ull, snp_t{"C"_dna5}, coverage_t{1, 1, 0, 0}});
    // jst.insert(event_t{4ull, snp_t{"T"_dna5}, coverage_t{0, 0, 0, 1}});

    std::vector matches = jstmap::search_queries_(pjst.bin_at(0), pattern_collection, 0.0);

    seqan3::debug_stream << "Report " << matches.size() << " matches:\n";
    for (auto const & match : matches)
        seqan3::debug_stream << "\t- match with: " << match.error_count << " errors, sequence = "
                             << match.jst_sequence << "\n";
    // EXPECT_EQ(results.size(), 10u);
}
