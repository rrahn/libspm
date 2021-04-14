// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <concepts>
#include <random>
#include <ranges>

#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/utility/detail/multi_invocable.hpp>

#include <jstmap/simulate/simulate_alignment.hpp>

#include <libjst/journaled_sequence_tree.hpp>

class jst_fuzzy_test : public ::testing::Test
{
public:

    using sequence_type = jstmap::sequence_t;
    using alignment_type = jstmap::alignment_t;
    using jst_type = libjst::journaled_sequence_tree<sequence_type>;

    sequence_type reference{};
    std::vector<alignment_type> alignments{};

    void SetUp() override
    {
        // ----------------------------------------------------------------------------
        // Initialise the random device.
        // ----------------------------------------------------------------------------

        std::random_device rd{};  //Will be used to obtain a seed for the random number engine
        uint64_t const seed = rd(); // We want to generate one seed for the entire test suite and only once per execution.
        std::mt19937_64 random_engine{seed};

        // ----------------------------------------------------------------------------
        // Generate the reference sequence.
        // ----------------------------------------------------------------------------

        seqan3::test::random_sequence_generator<sequence_type> sequence_generator{1000, 200};
        reference = sequence_generator(random_engine);

        // ----------------------------------------------------------------------------
        // Generate the alignments.
        // ----------------------------------------------------------------------------

        std::uniform_int_distribution<size_t> sequence_count_distribution(0, 100);
        std::uniform_real_distribution<double> error_rate_distribution(0.0, 0.5);

        double const sequence_count = sequence_count_distribution(random_engine);
        double const error_rate = error_rate_distribution(random_engine);

        alignments.resize(sequence_count);
        std::ranges::generate(alignments, [&] () { return jstmap::simulate_alignment(reference, error_rate); });

        // ----------------------------------------------------------------------------
        // Print simulation parameters
        // ----------------------------------------------------------------------------

        std::cout << "Simulation parameter:\n"
                  << "\t- Seed: " << seed << "\n"
                  << "\t- Reference size: " << reference.size() << "\n"
                  << "\t- Sequence count: " << sequence_count << "\n"
                  << "\t- Error rate: " << error_rate << "\n\n";
    }

    auto target_sequence(size_t const idx)
    {
        // Strip all gaps from the target sequence.
        return alignments[idx].second | std::views::filter([] (auto const & gapped_alphabet)
        {
            return gapped_alphabet != seqan3::gap{};
        });
    }
};

TEST_F(jst_fuzzy_test, jst_construction)
{
    jst_type jst{sequence_type{reference}};

    std::ranges::for_each(alignments, [&] (auto && alignment) { jst.add(alignment); });

    for (unsigned idx = 0; idx < alignments.size(); ++idx)
        EXPECT_RANGE_EQ(jst.sequence_at(idx), target_sequence(idx));
}
