// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <stack>
#include <string>



#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/stats.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/rcms/dna_compressed_multisequence.hpp>
#include <libjst/rcms/rcs_store.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

#include "../mock/rcs_store_mock.hpp"

namespace jst::test::sequence_tree_stats {

using source_t = std::string;
using variant_t = jst::test::variant<uint32_t, source_t, uint32_t, std::vector<uint32_t>>;

struct fixture {
    source_t source{};
    std::vector<variant_t> variants{};
    uint32_t coverage_size{};
    uint32_t window_size{};
    libjst::tree_stats expected_stats{};

    template <typename stream_t, typename this_t>
        requires std::same_as<std::remove_cvref_t<this_t>, fixture>
    friend stream_t & operator<<(stream_t & stream, this_t &&) {
        stream << "fixture";
        return stream;
    }
};

struct test : public ::testing::TestWithParam<fixture> {
    using coverage_type = libjst::bit_coverage<uint32_t>;
    using coverage_domain_type = libjst::coverage_domain_t<coverage_type>;

    using cms_t = libjst::dna_compressed_multisequence<source_t, coverage_type>;
    using cms_value_t = std::ranges::range_value_t<cms_t>;
    using rcs_store_t = libjst::rcs_store<source_t, cms_t>;
    rcs_store_t _mock;

    void SetUp() override {
        _mock = rcs_store_t{GetParam().source, GetParam().coverage_size};
        coverage_domain_type domain = _mock.variants().coverage_domain();

        std::ranges::for_each(GetParam().variants, [&] (auto var) {
            _mock.add(cms_value_t{libjst::breakpoint{var.position, var.deletion},
                                  var.insertion,
                                  coverage_type{var.coverage, domain}});
        });
    }

    rcs_store_t const & get_mock() const noexcept{
        return _mock;
    }
};

} // namespace jst::test::sequence_tree_stats

using namespace std::literals;

using fixture = jst::test::sequence_tree_stats::fixture;
using variant_t = jst::test::sequence_tree_stats::variant_t;
struct sequence_tree_stats : public jst::test::sequence_tree_stats::test
{
    using typename jst::test::sequence_tree_stats::test::rcs_store_t;
    using jst::test::sequence_tree_stats::test::get_mock;
    using jst::test::sequence_tree_stats::test::GetParam;

    auto make_tree() const noexcept {
        auto const & rcs_mock = get_mock();

        return libjst::make_volatile(rcs_mock) | libjst::labelled()
                                               | libjst::coloured()
                                               | libjst::trim(GetParam().window_size)
                                               | libjst::prune()
                                               | libjst::merge();
    }
};

// ----------------------------------------------------------------------------
// Test case definitions
// ----------------------------------------------------------------------------

TEST_P(sequence_tree_stats, node_count) {
    auto actual_stats = libjst::stats(make_tree());
    EXPECT_EQ(actual_stats.node_count, GetParam().expected_stats.node_count);
}

TEST_P(sequence_tree_stats, subtree_count) {
    auto actual_stats = libjst::stats(make_tree());
    EXPECT_EQ(actual_stats.subtree_count, GetParam().expected_stats.subtree_count);
}

TEST_P(sequence_tree_stats, leaf_count) {
    auto actual_stats = libjst::stats(make_tree());
    EXPECT_EQ(actual_stats.leaf_count, GetParam().expected_stats.leaf_count);
}

TEST_P(sequence_tree_stats, symbol_count) {
    auto actual_stats = libjst::stats(make_tree());
    EXPECT_EQ(actual_stats.symbol_count, GetParam().expected_stats.symbol_count);
}

TEST_P(sequence_tree_stats, max_subtree_depth) {
    auto actual_stats = libjst::stats(make_tree());
    EXPECT_EQ(actual_stats.max_subtree_depth, GetParam().expected_stats.max_subtree_depth);
}

TEST_P(sequence_tree_stats, subtree_depths) {
    auto actual_stats = libjst::stats(make_tree());
    EXPECT_TRUE(std::ranges::equal(actual_stats.subtree_depths, GetParam().expected_stats.subtree_depths));
}

// ----------------------------------------------------------------------------
// Test values
// ----------------------------------------------------------------------------

using namespace std::literals;

INSTANTIATE_TEST_SUITE_P(no_variants, sequence_tree_stats, testing::Values(fixture{
    .source{"AAAAGGGG"s},
    .variants{},
    .coverage_size{4},
    .window_size{4},
    .expected_stats{.node_count = 1,
                    .subtree_count = 0,
                    .leaf_count = 1,
                    .symbol_count = 8,
                    .max_subtree_depth = 0,
                    .subtree_depths{}}
}));

INSTANTIATE_TEST_SUITE_P(single_variant_first_base, sequence_tree_stats, testing::Values(fixture{
    .source{"AAAAGGGG"s},
    .variants{
        variant_t{.position{0}, .insertion{"C"s}, .deletion{1}, .coverage{0}}
    },
    .coverage_size{4},
    .window_size{4},
    .expected_stats{.node_count = 3,
                    .subtree_count = 1,
                    .leaf_count = 2,
                    .symbol_count = 8 + 5,
                    .max_subtree_depth = 1,
                    .subtree_depths{1}}
}));

INSTANTIATE_TEST_SUITE_P(single_variant_last_base, sequence_tree_stats, testing::Values(fixture{
    .source{"AAAAGGGG"s},
    .variants{
        variant_t{.position{7}, .insertion{"C"s}, .deletion{1}, .coverage{0}}
    },
    .coverage_size{4},
    .window_size{4},
    .expected_stats{.node_count = 3,
                    .subtree_count = 1,
                    .leaf_count = 2,
                    .symbol_count = 8 + 1,
                    .max_subtree_depth = 1,
                    .subtree_depths{1}}
}));

INSTANTIATE_TEST_SUITE_P(single_variant_middle, sequence_tree_stats, testing::Values(fixture{
    .source{"AAAAGGGG"s},
    .variants{
        variant_t{.position{4}, .insertion{"C"s}, .deletion{1}, .coverage{0}}
    },
    .coverage_size{4},
    .window_size{4},
    .expected_stats{.node_count = 3,
                    .subtree_count = 1,
                    .leaf_count = 2,
                    .symbol_count = 8 + 4,
                    .max_subtree_depth = 1,
                    .subtree_depths{1}}
}));

INSTANTIATE_TEST_SUITE_P(two_variants_non_overlapping, sequence_tree_stats, testing::Values(fixture{
         //  01234567
    .source{"AAAAGGGG"s},
    .variants{variant_t{.position{1}, .insertion{"C"s}, .deletion{1}, .coverage{0,1}},
              variant_t{.position{6}, .insertion{"T"s}, .deletion{1}, .coverage{0,2}}},
    .coverage_size{4},
    .window_size{4},
    .expected_stats{.node_count = 5,
                    .subtree_count = 2,
                    .leaf_count = 3,
                    .symbol_count = 8 + 5 + 2,
                    .max_subtree_depth = 1,
                    .subtree_depths{1, 1}}
}));

INSTANTIATE_TEST_SUITE_P(two_variants_overlapping, sequence_tree_stats, testing::Values(fixture{
         //  01234567
    .source{"AAAAGGGG"s},
    .variants{variant_t{.position{1}, .insertion{"C"s}, .deletion{1}, .coverage{0, 1}},
              variant_t{.position{4}, .insertion{"T"s}, .deletion{1}, .coverage{0, 2}}},
    .coverage_size{4},
    .window_size{4},
    .expected_stats{.node_count = 7,
                    .subtree_count = 2,
                    .leaf_count = 4,
                    .symbol_count = 8 + (5 + 2) + 4,
                    .max_subtree_depth = 2,
                    .subtree_depths{2, 1}}
}));

INSTANTIATE_TEST_SUITE_P(two_variants_overlapping_same_position, sequence_tree_stats, testing::Values(fixture{
         //  01234567
    .source{"AAAAGGGG"s},
    .variants{variant_t{.position{3}, .insertion{"C"s}, .deletion{1}, .coverage{0, 1}},
              variant_t{.position{3}, .insertion{"T"s}, .deletion{1}, .coverage{0, 2}}},
    .coverage_size{4},
    .window_size{4},
    .expected_stats{.node_count = 5,
                    .subtree_count = 2,
                    .leaf_count = 3,
                    .symbol_count = 8 + 5 + 5,
                    .max_subtree_depth = 1,
                    .subtree_depths{1, 1}}
}));
