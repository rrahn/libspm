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

#include <libcontrib/seqan/alphabet.hpp>

#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/partial_tree.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/rcms/compressed_multisequence.hpp>
#include <libjst/rcms/rcs_store.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

#include "../mock/rcs_store_mock.hpp"

namespace jst::test::partial_sequence_tree {

using source_t = std::vector<jst::contrib::dna4>;
using variant_t = jst::test::variant<uint32_t, source_t, uint32_t, std::vector<uint32_t>>;

struct fixture {
    source_t source{};
    std::vector<variant_t> variants{};
    uint32_t coverage_size{};
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
    using coverage_type = libjst::bit_coverage<uint32_t>;
    using coverage_domain_type = libjst::coverage_domain_t<coverage_type>;

    using cms_t = libjst::compressed_multisequence<source_t, coverage_type>;
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
        libjst::partial_tree partial_mock{rcs_mock, GetParam().bin_offset, GetParam().bin_size};

        return std::move(partial_mock) | libjst::labelled()
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

    using tree_t = decltype(tree);

    using node_t = libjst::tree_node_t<tree_t>;
    using cargo_t = libjst::tree_label_t<tree_t>;

    auto to_string = [] (auto seq) -> std::string {
        std::string str;
        for (char c : seq)
            str.push_back(c);
        return str;
    };

    std::vector<std::string> actual_labels{};
    std::stack<node_t> path{};
    path.push(libjst::root(tree));

    std::cout << "Labels: ";
    while (!path.empty()) {
        node_t p = std::move(path.top());
        path.pop();
        cargo_t label = *p;
        if (!std::ranges::empty(label.sequence())) {
            actual_labels.push_back(to_string(label.sequence()));
            std::cout << actual_labels.back() << " " << std::flush;
        }

        if (auto c_ref = p.next_ref(); c_ref.has_value()) {
            path.push(std::move(*c_ref));
        }
        if (auto c_alt = p.next_alt(); c_alt.has_value()) {
            path.push(std::move(*c_alt));
        }
    }
    std::cout << "\n";

    std::ptrdiff_t expected_count = std::ranges::ssize(GetParam().expected_labels);
    std::ptrdiff_t actual_count = std::ranges::ssize(actual_labels);
    EXPECT_EQ(expected_count, actual_count);
    for (std::ptrdiff_t i = 0; i < std::min(expected_count, actual_count); ++i)
        EXPECT_EQ(to_string(GetParam().expected_labels[i]), actual_labels[i]) << i;
}

// ----------------------------------------------------------------------------
// Test values
// ----------------------------------------------------------------------------
using jst::contrib::operator""_dna4;

INSTANTIATE_TEST_SUITE_P(no_variant_unbound, partial_sequence_tree_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{},
    .coverage_size{4},
    .bin_offset{0},
    .bin_size{8},
    .window_size{4},
    .expected_labels{"AAAAGGGG"_dna4}
}));

INSTANTIATE_TEST_SUITE_P(no_variant_left_bound, partial_sequence_tree_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{},
    .coverage_size{4},
    .bin_offset{2},
    .bin_size{6},
    .window_size{4},
    .expected_labels{"AAGGGG"_dna4}
}));

INSTANTIATE_TEST_SUITE_P(no_variant_right_bound, partial_sequence_tree_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{},
    .coverage_size{4},
    .bin_offset{0},
    .bin_size{6},
    .window_size{4},
    .expected_labels{"AAAAGG"_dna4, "GG"_dna4}
}));

INSTANTIATE_TEST_SUITE_P(no_variant_left_and_right_bound, partial_sequence_tree_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{},
    .coverage_size{4},
    .bin_offset{2},
    .bin_size{4},
    .window_size{4},
    .expected_labels{"AAGG"_dna4, "GG"_dna4}
}));

INSTANTIATE_TEST_SUITE_P(no_variant_left_and_right_bound_single, partial_sequence_tree_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{},
    .coverage_size{4},
    .bin_offset{2},
    .bin_size{1},
    .window_size{4},
    .expected_labels{"A"_dna4, "AGGG"_dna4}
}));

INSTANTIATE_TEST_SUITE_P(two_variants_unbound, partial_sequence_tree_test, testing::Values(fixture{
         //  01234567
    .source{"AAAAGGGG"_dna4},
    .variants{variant_t{.position{1}, .insertion{"C"_dna4}, .deletion{1}, .coverage{0,1}},
              variant_t{.position{5}, .insertion{"T"_dna4}, .deletion{1}, .coverage{0,2}}},
    .coverage_size{4},
    .bin_offset{0},
    .bin_size{8},
    .window_size{4},
    .expected_labels{"A"_dna4,
                       "CAAG"_dna4,
                            "T"_dna4,
                            "G"_dna4,
                       "AAAG"_dna4,
                            "TGG"_dna4,
                            "GGG"_dna4,}
}));

INSTANTIATE_TEST_SUITE_P(two_variants_left_bound, partial_sequence_tree_test, testing::Values(fixture{
         //  01234567
    .source{"AAAAGGGG"_dna4},
    .variants{variant_t{.position{1}, .insertion{"C"_dna4}, .deletion{1}, .coverage{0,1}},
              variant_t{.position{5}, .insertion{"T"_dna4}, .deletion{1}, .coverage{0,2}}},
    .coverage_size{4},
    .bin_offset{1},
    .bin_size{7},
    .window_size{4},
    .expected_labels{"CAAG"_dna4,
                          "T"_dna4,
                          "G"_dna4,
                     "AAAG"_dna4,
                          "TGG"_dna4,
                          "GGG"_dna4,}
}));

INSTANTIATE_TEST_SUITE_P(two_variants_right_bound, partial_sequence_tree_test, testing::Values(fixture{
         //  01234567
    .source{"AAAAGGGG"_dna4},
    .variants{variant_t{.position{1}, .insertion{"C"_dna4}, .deletion{1}, .coverage{0,1}},
              variant_t{.position{5}, .insertion{"T"_dna4}, .deletion{1}, .coverage{0,2}}},
    .coverage_size{4},
    .bin_offset{0},
    .bin_size{5},
    .window_size{4},
    .expected_labels{"A"_dna4,
                      "CAAG"_dna4,
                          "T"_dna4,
                          "G"_dna4,
                      "AAAG"_dna4, // end after here, vvvv overlaps
                          "TGG"_dna4,
                          "GGG"_dna4}
}));

INSTANTIATE_TEST_SUITE_P(two_variants_left_and_right_bound_inclusive, partial_sequence_tree_test, testing::Values(fixture{
         //  01234567
    .source{"AAAAGGGG"_dna4},
    .variants{variant_t{.position{1}, .insertion{"C"_dna4}, .deletion{1}, .coverage{0,1}},
              variant_t{.position{4}, .insertion{"T"_dna4}, .deletion{1}, .coverage{0,2}}},
    .coverage_size{4},
    .bin_offset{1},
    .bin_size{4},
    .window_size{4},
    .expected_labels{"CAA"_dna4,
                        "TG"_dna4,
                        "GG"_dna4,
                     "AAA"_dna4,
                        "TGGG"_dna4,
                        "G"_dna4,
                         "GGG"_dna4}
}));

INSTANTIATE_TEST_SUITE_P(two_variants_left_and_right_bound_exclusive, partial_sequence_tree_test, testing::Values(fixture{
         //  01234567
    .source{"AAAAGGGG"_dna4},
    .variants{variant_t{.position{1}, .insertion{"C"_dna4}, .deletion{1}, .coverage{0,1}},
              variant_t{.position{5}, .insertion{"T"_dna4}, .deletion{1}, .coverage{0,2}}},
    .coverage_size{4},
    .bin_offset{2},
    .bin_size{3},
    .window_size{4},
    .expected_labels{"AAG"_dna4,
                        "TGG"_dna4,
                        "GGG"_dna4}
}));

