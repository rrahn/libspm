// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <concepts>
#include <algorithm>
#include <stack>
#include <string>

#include <seqan3/test/expect_range_eq.hpp>

#include <libcontrib/seqan/alphabet.hpp>

#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/rcms/compressed_multisequence.hpp>
#include <libjst/rcms/rcs_store.hpp>
#include <libjst/rcms/rcs_store_reversed.hpp>

#include "../mock/rcs_store_mock.hpp"

namespace jst::test::merged_tree_reverse {

using source_t = std::vector<jst::contrib::dna4>;
using variant_t = jst::test::variant<uint32_t, source_t, uint32_t, std::vector<uint32_t>>;

struct fixture {
    source_t source{};
    uint32_t coverage_size{4};
    std::vector<variant_t> variants{};
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
    using rcs_store_reverse_t = libjst::rcs_store_reversed<cms_t>;

    rcs_store_t _mock;
    std::unique_ptr<rcs_store_reverse_t> _reversed_mock{};

    void SetUp() override {
        _mock = rcs_store_t{GetParam().source, GetParam().coverage_size};
        coverage_domain_type domain = _mock.variants().coverage_domain();

        std::ranges::for_each(GetParam().variants, [&] (auto var) {
            _mock.add(cms_value_t{libjst::breakpoint{var.position, var.deletion},
                                  var.insertion,
                                  coverage_type{var.coverage, domain}});
        });
        _reversed_mock = std::make_unique<rcs_store_reverse_t>(_mock.variants());
    }

    rcs_store_reverse_t const & get_mock() const noexcept{
        return *_reversed_mock;
    }
};

} // namespace jst::test::merged_tree_reverse

using namespace std::literals;

using fixture = jst::test::merged_tree_reverse::fixture;
using variant_t = jst::test::merged_tree_reverse::variant_t;
struct merged_tree_reverse_test : public jst::test::merged_tree_reverse::test
{
    using typename jst::test::merged_tree_reverse::test::rcs_store_t;
    using jst::test::merged_tree_reverse::test::get_mock;
    using jst::test::merged_tree_reverse::test::GetParam;

    auto make_tree() const noexcept {
        auto const & rcs_mock = get_mock();
        return libjst::volatile_tree{rcs_mock} | libjst::labelled() | libjst::merge();
    }
};

// ----------------------------------------------------------------------------
// Test case definitions
// ----------------------------------------------------------------------------

TEST_P(merged_tree_reverse_test, root_sink) {
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

    node_t r = libjst::root(tree);

    std::vector<std::string> actual_labels{};
    std::stack<node_t> path{};
    path.push(r);

    std::cout << "Labels: ";
    while (!path.empty()) {
        node_t p = std::move(path.top());
        path.pop();
        cargo_t label = *p;
        std::cout << to_string(label.sequence()) << " " << std::flush;
        actual_labels.push_back(to_string(label.sequence()));

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

INSTANTIATE_TEST_SUITE_P(no_variant, merged_tree_reverse_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{},
    .expected_labels{"GGGGAAAA"_dna4}
}));

INSTANTIATE_TEST_SUITE_P(snv0, merged_tree_reverse_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{
        variant_t{.position{0}, .insertion{"C"_dna4}, .deletion{1}, .coverage{0}}
    },
    .expected_labels{"GGGGAAA"_dna4, "C"_dna4, "A"_dna4}
}));

INSTANTIATE_TEST_SUITE_P(snv7, merged_tree_reverse_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{
        variant_t{.position{7}, .insertion{"C"_dna4}, .deletion{1}, .coverage{0}}
    },
    .expected_labels{""_dna4, "CGGGAAAA"_dna4, "GGGGAAAA"_dna4}
}));

INSTANTIATE_TEST_SUITE_P(snv4, merged_tree_reverse_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{
        variant_t{.position{4}, .insertion{"C"_dna4}, .deletion{1}, .coverage{0}}
    },
    .expected_labels{"GGG"_dna4, "CAAAA"_dna4, "GAAAA"_dna4}
}));

INSTANTIATE_TEST_SUITE_P(snv4_snv6, merged_tree_reverse_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{
        variant_t{.position{4}, .insertion{"C"_dna4}, .deletion{1}, .coverage{0}},
        variant_t{.position{6}, .insertion{"T"_dna4}, .deletion{1}, .coverage{0, 2}}
    },
    .expected_labels{"G"_dna4, "TG"_dna4, "CAAAA"_dna4,
                                          "GAAAA"_dna4,
                               "GG"_dna4, "CAAAA"_dna4,
                                          "GAAAA"_dna4}
}));

INSTANTIATE_TEST_SUITE_P(snv4_snv5, merged_tree_reverse_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{
        variant_t{.position{4}, .insertion{"C"_dna4}, .deletion{1}, .coverage{0}},
        variant_t{.position{5}, .insertion{"T"_dna4}, .deletion{1}, .coverage{0, 2}}
    },
    .expected_labels{"GG"_dna4, "T"_dna4, "CAAAA"_dna4,
                                          "GAAAA"_dna4,
                                "G"_dna4, "CAAAA"_dna4,
                                          "GAAAA"_dna4}
}));

INSTANTIATE_TEST_SUITE_P(snv4_snv4, merged_tree_reverse_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{
        variant_t{.position{4}, .insertion{"C"_dna4}, .deletion{1}, .coverage{0}},
        variant_t{.position{4}, .insertion{"T"_dna4}, .deletion{1}, .coverage{1, 2}}
    },
    .expected_labels{"GGG"_dna4, "TAAAA"_dna4,
                                  ""_dna4, "CAAAA"_dna4,
                                  "GAAAA"_dna4}
}));
