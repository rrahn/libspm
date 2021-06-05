// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>
#include <set>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include <libjst/journaled_sequence_tree.hpp>
#include <libjst/search/pigeonhole_filter.hpp>
#include <libjst/search/state_manager_stack.hpp>

#include "../test_utility.hpp" // load_jst

using namespace std::literals;

using sequence_t = seqan3::dna5_vector;

struct pigeonhole_filter_test : public ::testing::Test
{
};

TEST(pigeonhole_filter_test, construction)
{
    // using seqan3::operator""_dna5;

    // using jst_t = libjst::journaled_sequence_tree<sequence_t>;

    // seqan::StringSet<sequence_t> pattern_collection{};
    // seqan::appendValue(pattern_collection, "acgtaacgtaacgtagacga"_dna5);
    // seqan::appendValue(pattern_collection, "acgtacgactacgtacgact"_dna5);
    // seqan::appendValue(pattern_collection, "acgtacgactagcgactacg"_dna5);

    // std::filesystem::path jst_file{DATADIR"sim_refx5.jst"};
    // jst_t jst = libjst::test::load_jst<jst_t>(jst_file);

    // using state_t = typename decltype(libjst::pigeonhole_filter{pattern_collection})::state_type;
    // // filter set error rate.
    // // set state manager from outside.
    // libjst::pigeonhole_filter filter{pattern_collection, 0.0, libjst::search_state_manager_stack<state_t>{}};
    // auto jst_range_agent = jst.range_agent(filter.qgram_size(), filter.state_manager()); // already pushing a branch.
    // // sequence_t haystack{"acgtaacgtaacgtagacgaacgtacgactacgtacgactacgtacgactagcgactacg"_dna5};
    // filter(jst_range_agent, [&] (auto const & hit, auto const & haystack_it)
    // {
    //     std::cout << "hit: " << hit << " at " << haystack_it.coordinate() << "\n";
    // });
}
