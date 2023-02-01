// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <string>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/views/slice.hpp>

#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/partial_tree.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

#include "../mock/rcs_store_mock.hpp"

namespace jst::test::partial_sequence_tree {

using source_t = std::string;
using variant_t = jst::test::variant<libjst::breakpoint, source_t, int, libjst::bit_vector<>>;

struct expected_root {
    source_t sequence{};
};

struct fixture {
    source_t source{};
    std::vector<variant_t> variants{};
    std::size_t coverage_size{};
    uint32_t bin_offset{};
    uint32_t bin_size{};
    uint32_t window_size{};
    std::vector<source_t> expected_labels{};

    template <typename stream_t, typename this_t>
        requires std::same_as<std::remove_cvref_t<this_t>, fixture>
    friend stream_t & operator<<(stream_t & stream, this_t &&) {
        stream << "fixture";
        return stream;
    }
};

struct test : public ::testing::TestWithParam<fixture> {

    using rcs_store_t = mock_store<source_t>;

    rcs_store_t _mock;

    void SetUp() override {
        _mock = rcs_store_t{GetParam().source, GetParam().coverage_size};
        std::ranges::for_each(GetParam().variants, [&] (auto var) {
            assert(std::ranges::size(libjst::coverage(var)) == _mock.size());
            _mock.insert(std::move(var));
        });
    }

    rcs_store_t const & get_mock() const noexcept{
        return _mock;
    }
};

} // namespace jst::test::partial_sequence_tree

using namespace std::literals;

using fixture = jst::test::partial_sequence_tree::fixture;
using variant_t = jst::test::partial_sequence_tree::variant_t;
struct partial_sequence_tree_test : public jst::test::partial_sequence_tree::test
{
    using typename jst::test::partial_sequence_tree::test::rcs_store_t;
    using jst::test::partial_sequence_tree::test::get_mock;
    using jst::test::partial_sequence_tree::test::GetParam;

    auto make_tree() const noexcept {
        auto const & rcs_mock = get_mock();
        libjst::volatile_tree mock_tree{rcs_mock};
        libjst::partial_tree partial_mock{std::move(mock_tree), GetParam().bin_offset, GetParam().bin_size};

        return std::move(partial_mock) | libjst::labelled<libjst::sequence_label_kind::root_path>()
                                       | libjst::coloured()
                                       | libjst::trim(GetParam().window_size)
                                       | libjst::merge();
    }
};

// ----------------------------------------------------------------------------
// Test case definitions
// ----------------------------------------------------------------------------

TEST_P(partial_sequence_tree_test, traverse) {
    auto tree = make_tree();

    libjst::tree_traverser_base path{tree};
    auto expected_it = std::ranges::begin(GetParam().expected_labels);
    for (auto it = path.begin(); it != path.end(); ++it) {
        auto && label = *it;
        if (std::ranges::empty(label.sequence()))
            continue;

        EXPECT_TRUE(expected_it != std::ranges::end(GetParam().expected_labels));
        EXPECT_RANGE_EQ(*expected_it, label.sequence());
        if (!std::ranges::equal(*expected_it, label.sequence()))
            std::cout << "Label pos: " << std::ranges::distance(std::ranges::begin(GetParam().expected_labels), expected_it) << "\n";
        ++expected_it;
    }
}

// ----------------------------------------------------------------------------
// Test values
// ----------------------------------------------------------------------------

INSTANTIATE_TEST_SUITE_P(no_variant_unbound, partial_sequence_tree_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{},
    .coverage_size{4},
    .bin_offset{0},
    .bin_size{8},
    .window_size{4},
    .expected_labels{"aaaabbbb"}
}));

INSTANTIATE_TEST_SUITE_P(no_variant_left_bound, partial_sequence_tree_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{},
    .coverage_size{4},
    .bin_offset{2},
    .bin_size{6},
    .window_size{4},
    .expected_labels{"aabbbb"}
}));

INSTANTIATE_TEST_SUITE_P(no_variant_right_bound, partial_sequence_tree_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{},
    .coverage_size{4},
    .bin_offset{0},
    .bin_size{6},
    .window_size{4},
    .expected_labels{"aaaabb", "bb"}
}));

INSTANTIATE_TEST_SUITE_P(no_variant_left_and_right_bound, partial_sequence_tree_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{},
    .coverage_size{4},
    .bin_offset{2},
    .bin_size{4},
    .window_size{4},
    .expected_labels{"aabb"}
}));

INSTANTIATE_TEST_SUITE_P(no_variant_left_and_right_bound_single, partial_sequence_tree_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{},
    .coverage_size{4},
    .bin_offset{2},
    .bin_size{1},
    .window_size{4},
    .expected_labels{"a"}
}));

INSTANTIATE_TEST_SUITE_P(two_variants_unbound, partial_sequence_tree_test, testing::Values(fixture{
         //  01234567
    .source{"aaaabbbb"},
    .variants{variant_t{.position{1}, .insertion{"I"}, .deletion{1}, .coverage{1,1,0,0}},
              variant_t{.position{5}, .insertion{"J"}, .deletion{1}, .coverage{1,0,1,0}}},
    .coverage_size{4},
    .bin_offset{0},
    .bin_size{8},
    .window_size{4},
    .expected_labels{"a",
                       "Iaab",
                            "J",
                            "b",
                       "aaab",
                            "Jbb",
                            "bbb",}
}));


INSTANTIATE_TEST_SUITE_P(two_variants_left_bound, partial_sequence_tree_test, testing::Values(fixture{
         //  01234567
    .source{"aaaabbbb"},
    .variants{variant_t{.position{1}, .insertion{"I"}, .deletion{1}, .coverage{1,1,0,0}},
              variant_t{.position{5}, .insertion{"J"}, .deletion{1}, .coverage{1,0,1,0}}},
    .coverage_size{4},
    .bin_offset{1},
    .bin_size{7},
    .window_size{4},
    .expected_labels{"Iaab",
                          "J",
                          "b",
                     "aaab",
                          "Jbb",
                          "bbb",}
}));

INSTANTIATE_TEST_SUITE_P(two_variants_right_bound, partial_sequence_tree_test, testing::Values(fixture{
         //  01234567
    .source{"aaaabbbb"},
    .variants{variant_t{.position{1}, .insertion{"I"}, .deletion{1}, .coverage{1,1,0,0}},
              variant_t{.position{5}, .insertion{"J"}, .deletion{1}, .coverage{1,0,1,0}}},
    .coverage_size{4},
    .bin_offset{0},
    .bin_size{5},
    .window_size{4},
    .expected_labels{"a",
                       "Iaab",
                            "J",
                            "b",
                       "aaab"}
}));

INSTANTIATE_TEST_SUITE_P(two_variants_left_and_right_bound_inclusive, partial_sequence_tree_test, testing::Values(fixture{
         //  01234567
    .source{"aaaabbbb"},
    .variants{variant_t{.position{1}, .insertion{"I"}, .deletion{1}, .coverage{1,1,0,0}},
              variant_t{.position{4}, .insertion{"J"}, .deletion{1}, .coverage{1,0,1,0}}},
    .coverage_size{4},
    .bin_offset{1},
    .bin_size{4},
    .window_size{4},
    .expected_labels{"Iaa",
                        "Jb",
                        "bb",
                     "aaa",
                        "Jbbb",
                        "b"}
}));

INSTANTIATE_TEST_SUITE_P(two_variants_left_and_right_bound_exclusive, partial_sequence_tree_test, testing::Values(fixture{
         //  01234567
    .source{"aaaabbbb"},
    .variants{variant_t{.position{1}, .insertion{"I"}, .deletion{1}, .coverage{1,1,0,0}},
              variant_t{.position{5}, .insertion{"J"}, .deletion{1}, .coverage{1,0,1,0}}},
    .coverage_size{4},
    .bin_offset{2},
    .bin_size{3},
    .window_size{4},
    .expected_labels{"aab"}
}));

