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
#include <seqan3/utility/views/slice.hpp>

#include <libcontrib/seqan/alphabet.hpp>

#include <libjst/sequence_tree/empty_label.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/rcms/compressed_multisequence.hpp>
#include <libjst/rcms/rcs_store.hpp>
#include <libjst/rcms/rcs_store_reversed.hpp>

#include "../mock/rcs_store_mock.hpp"

namespace jst::test::volatile_tree_reverse {

using source_t = std::vector<jst::contrib::dna4>;
using variant_t = jst::test::variant<uint32_t, source_t, uint32_t, std::vector<uint32_t>>;

struct fixture {
    source_t source{};
    uint32_t coverage_size{4};
    std::vector<variant_t> variants{};
    std::vector<uint32_t> expected_traversal{};

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
    using reversed_rcs_store_t = libjst::rcs_store_reversed<cms_t>;

    rcs_store_t _mock;
    std::unique_ptr<reversed_rcs_store_t> _reversed_mock{};

    void SetUp() override {
        _mock = rcs_store_t{GetParam().source, GetParam().coverage_size};
        coverage_domain_type domain = _mock.variants().coverage_domain();

        std::ranges::for_each(GetParam().variants, [&] (auto var) {
            _mock.add(cms_value_t{libjst::breakpoint{var.position, var.deletion},
                                  var.insertion,
                                  coverage_type{var.coverage, domain}});
        });
        _reversed_mock = std::make_unique<reversed_rcs_store_t>(_mock.variants());
    }

    reversed_rcs_store_t const & get_mock() const noexcept{
        return *_reversed_mock;
    }
};

} // namespace jst::test::volatile_tree_reverse

using namespace std::literals;

using fixture = jst::test::volatile_tree_reverse::fixture;
using variant_t = jst::test::volatile_tree_reverse::variant_t;
struct volatile_tree_reverse_test : public jst::test::volatile_tree_reverse::test
{
    using typename jst::test::volatile_tree_reverse::test::rcs_store_t;
    using jst::test::volatile_tree_reverse::test::get_mock;
    using jst::test::volatile_tree_reverse::test::GetParam;

    auto make_tree() const noexcept {
        auto const & rcs_mock = get_mock();
        return libjst::volatile_tree{rcs_mock};
    }
};

// ----------------------------------------------------------------------------
// Test case definitions
// ----------------------------------------------------------------------------

TEST_P(volatile_tree_reverse_test, root_sink) {
    auto tree = make_tree();

    using tree_t = decltype(tree);

    using node_t = libjst::tree_node_t<tree_t>;
    using sink_t = libjst::tree_sink_t<tree_t>;
    using cargo_t = libjst::tree_label_t<tree_t>;

    EXPECT_TRUE((std::same_as<cargo_t, libjst::empty_label>));

    node_t r = libjst::root(tree);
    sink_t s = libjst::sink(tree);

    uint32_t visitor_id{0};
    std::vector<uint32_t> actual_traversal{};
    std::stack<node_t> path{};
    std::cout << "id: ";
    if (r != s) {
        std::cout << visitor_id << ", ";
        actual_traversal.push_back(visitor_id++);
        path.push(r);
    }
    while (!path.empty()) {
        node_t p = std::move(path.top());
        path.pop();

        if (auto c_ref = p.next_ref(); c_ref.has_value()) {
            std::cout << visitor_id << ", ";
            actual_traversal.push_back(visitor_id++);
            path.push(std::move(*c_ref));
        }
        if (auto c_alt = p.next_alt(); c_alt.has_value()) {
            std::cout << visitor_id << ", ";
            actual_traversal.push_back(visitor_id++);
            path.push(std::move(*c_alt));
        }
    }
    std::cout << "\n";

    EXPECT_RANGE_EQ(GetParam().expected_traversal, actual_traversal);
}

// ----------------------------------------------------------------------------
// Test values
// ----------------------------------------------------------------------------
using jst::contrib::operator""_dna4;

INSTANTIATE_TEST_SUITE_P(no_variant, volatile_tree_reverse_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{},
    .expected_traversal{}
}));

INSTANTIATE_TEST_SUITE_P(snv0, volatile_tree_reverse_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{
        variant_t{.position{0}, .insertion{"C"_dna4}, .deletion{1}, .coverage{0}}
    },
    .expected_traversal{0, 1, 2, 3}
}));

INSTANTIATE_TEST_SUITE_P(snv7, volatile_tree_reverse_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{
        variant_t{.position{7}, .insertion{"C"_dna4}, .deletion{1}, .coverage{0}}
    },
    .expected_traversal{0, 1, 2, 3}
}));

INSTANTIATE_TEST_SUITE_P(snv4, volatile_tree_reverse_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{
        variant_t{.position{4}, .insertion{"C"_dna4}, .deletion{1}, .coverage{0}}
    },
    .expected_traversal{0, 1, 2, 3}
}));

INSTANTIATE_TEST_SUITE_P(snv4_snv6, volatile_tree_reverse_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{
        variant_t{.position{4}, .insertion{"C"_dna4}, .deletion{1}, .coverage{0}},
        variant_t{.position{6}, .insertion{"T"_dna4}, .deletion{1}, .coverage{0, 2}}
    },
    .expected_traversal{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
}));

INSTANTIATE_TEST_SUITE_P(snv4_snv5, volatile_tree_reverse_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{
        variant_t{.position{4}, .insertion{"C"_dna4}, .deletion{1}, .coverage{0}},
        variant_t{.position{5}, .insertion{"T"_dna4}, .deletion{1}, .coverage{0, 2}}
    },
    .expected_traversal{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
}));

INSTANTIATE_TEST_SUITE_P(snv4_snv4, volatile_tree_reverse_test, testing::Values(fixture{
    .source{"AAAAGGGG"_dna4},
    .variants{
        variant_t{.position{4}, .insertion{"C"_dna4}, .deletion{1}, .coverage{0}},
        variant_t{.position{4}, .insertion{"T"_dna4}, .deletion{1}, .coverage{1, 2}}
    },
    .expected_traversal{0, 1, 2, 3, 4, 5, 6}
}));
