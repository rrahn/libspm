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
#include <string>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/range/views/char_to.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/utility/detail/multi_invocable.hpp>

#include <jstmap/simulate/simulate_alignment.hpp>

#include <libjst/set.hpp>

#include "test_utility.hpp"
class jst_fuzzy_test : public ::testing::Test,
                       public libjst::test::jst_context_map_fixture
{
public:

    using sequence_type = jstmap::raw_sequence_t;
    using alphabet_type = std::ranges::range_value_t<sequence_type>;
    using alignment_type = jstmap::alignment_t;
    using jst_type = libjst::journaled_sequence_tree<sequence_type>;

    sequence_type reference{};
    std::vector<alignment_type> alignments{};
    std::vector<std::string> target_sequences{};

    inline static uint64_t seed{};
    inline static size_t reference_size{};
    inline static size_t context_size{};
    inline static size_t sequence_count{};
    inline static double error_rate{};

    static void SetUpTestCase()
    {
        // ----------------------------------------------------------------------------
        // Initialise the seed globally
        // ----------------------------------------------------------------------------

        static std::random_device rd{};
        seed = rd();
        std::cout << "Seed: " << seed << "\n";
    }

    void SetUp() override
    {
        // ----------------------------------------------------------------------------
        // Generate the reference sequence and context.
        // ----------------------------------------------------------------------------

        std::mt19937_64 random_engine{seed};

        seqan3::test::random_sequence_generator<sequence_type> sequence_generator{800, 200};
        reference = sequence_generator(random_engine);
        reference_size = reference.size();

        std::uniform_int_distribution<size_t> context_size_distribution(1, 50);
        context_size = context_size_distribution(random_engine);

        // ----------------------------------------------------------------------------
        // Generate the alignments.
        // ----------------------------------------------------------------------------

        std::uniform_int_distribution<size_t> sequence_count_distribution(0, 50);
        std::uniform_real_distribution<double> error_rate_distribution(0.0, 0.1);

        sequence_count = sequence_count_distribution(random_engine);
        error_rate = error_rate_distribution(random_engine);

        alignments.resize(sequence_count);
        std::ranges::generate(alignments, [&] () { return jstmap::simulate_alignment(reference, error_rate); });

        // ----------------------------------------------------------------------------
        // Extract the target sequences and generate context map.
        // ----------------------------------------------------------------------------

        std::ranges::for_each(alignments, [&] (alignment_type const & alignment)
        {
            target_sequences.push_back(target_sequence(alignment));
        });

        generate_context_map(context_size, target_sequences);

    }

    static void TearDownTestCase()
    {
        // ----------------------------------------------------------------------------
        // Print simulation parameters
        // ----------------------------------------------------------------------------

        std::cout << "Simulation parameter:\n"
                  << "\t- Seed: " << seed << "\n"
                  << "\t- Reference size: " << reference_size << "\n"
                  << "\t- Context size: " << context_size << "\n"
                  << "\t- Sequence count: " << sequence_count << "\n"
                  << "\t- Error rate: " << error_rate << "\n\n";
    }

    auto target_sequence_at(size_t const idx) const
        -> decltype(std::declval<std::string const &>() | seqan3::views::char_to<alphabet_type>)
    {
        return target_sequences[idx] | seqan3::views::char_to<alphabet_type>;
    }

private:

    std::string target_sequence(alignment_type const & alignment) const
    {
        // Strip all gaps from the target sequence.
        return libjst::test::sequence_to_string(alignment.second | std::views::filter([] (auto const & gapped_alphabet)
        {
            return gapped_alphabet != seqan3::gap{};
        }));
    }
};

TEST_F(jst_fuzzy_test, jst_construction)
{
    jst_type jst{sequence_type{reference}};

    std::ranges::for_each(alignments, [&] (auto && alignment) { jst.add(alignment); });

    for (unsigned idx = 0; idx < jst.size(); ++idx)
        EXPECT_RANGE_EQ(jst.sequence_at(idx), target_sequence_at(idx));
}

TEST_F(jst_fuzzy_test, jst_context_enumeration)
{
    jst_type jst{sequence_type{reference}};
    std::ranges::for_each(alignments, [&] (auto && alignment) { jst.add(alignment); });

    auto context_enumerator = jst.context_enumerator(context_size);

    for (auto context_it = context_enumerator.begin(); context_it != context_enumerator.end(); ++context_it)
    {
        auto context = *context_it;
        std::string tmp = libjst::test::sequence_to_string(context);

        auto positions = jst.sequence_positions_at(context_it.coordinate());

        EXPECT_TRUE((this->context_positions_exist(tmp, positions))) << "context " << tmp;
    }

    // Verify that all unique contexts have been enumerated and that there is no unknown location.
    EXPECT_TRUE(this->all_contexts_enumerated());
    print_unvisited_contexts();

    EXPECT_TRUE(this->unknown_locations.empty());
    print_unknown_context_locations();
}
