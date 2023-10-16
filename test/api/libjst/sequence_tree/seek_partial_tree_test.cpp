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



#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/partial_tree.hpp>
#include <libjst/rcms/dna_compressed_multisequence.hpp>
#include <libjst/rcms/rcs_store.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

#include "../mock/rcs_store_mock.hpp"

namespace jst::test::seek_partial_tree {

using source_t = std::string;
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

} // namespace jst::test::seek_partial_tree

using namespace std::literals;

using fixture = jst::test::seek_partial_tree::fixture;
using variant_t = jst::test::seek_partial_tree::variant_t;
struct seek_partial_tree_test : public jst::test::seek_partial_tree::test
{
    using typename jst::test::seek_partial_tree::test::rcs_store_t;
    using jst::test::seek_partial_tree::test::get_mock;
    using jst::test::seek_partial_tree::test::GetParam;

    auto make_tree() const noexcept {
        auto const & rcs_mock = get_mock();
        return libjst::partial_tree{rcs_mock, GetParam().bin_offset, GetParam().bin_size}
                                               | libjst::labelled()
                                               | libjst::trim(GetParam().window_size)
                                               | libjst::merge();
    }
};

// ----------------------------------------------------------------------------
// Test case definitions
// ----------------------------------------------------------------------------

TEST_P(seek_partial_tree_test, seek) {
    auto tree = make_tree() | libjst::seek();

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
    // size_t counter{};
    while (!path.empty()) {
        node_t p = std::move(path.top());
        path.pop();
        cargo_t expected_label = *p;
        auto seek_pos = expected_label.position();
        // std::cout << counter++ << ": " << seek_pos << "\n";
        auto test_tree = make_tree() | libjst::seek();
        auto tmp = test_tree.seek(seek_pos);
        cargo_t actual_label = *tmp;
        EXPECT_EQ(to_string(actual_label.sequence()), to_string(expected_label.sequence()));
        std::cout << to_string(actual_label.sequence()) << " ";

        if (auto c_ref = p.next_ref(); c_ref.has_value()) {
            path.push(std::move(*c_ref));
        }
        if (auto c_alt = p.next_alt(); c_alt.has_value()) {
            path.push(std::move(*c_alt));
        }
    }
    std::cout << "\n";
}

// ----------------------------------------------------------------------------
// Test values
// ----------------------------------------------------------------------------

using namespace std::literals;

INSTANTIATE_TEST_SUITE_P(variant_on_partial_sink, seek_partial_tree_test, testing::Values(fixture{
    .source{"AAAACCCCGGGGTTTT"s},
    .variants{variant_t{.position{8}, .insertion{"A"s}, .deletion{1}, .coverage{0,1}}},
    .coverage_size{4},
    .bin_offset{4},
    .bin_size{4},
    .window_size{3},
    .expected_labels{"CCCC"s, "AGG"s, "GGG"s}
}));

INSTANTIATE_TEST_SUITE_P(variant_before_partial_sink, seek_partial_tree_test, testing::Values(fixture{
    .source{"AAAACCCCGGGGTTTT"s},
    .variants{variant_t{.position{7}, .insertion{"A"s}, .deletion{1}, .coverage{0,1}}},
    .coverage_size{4},
    .bin_offset{4},
    .bin_size{4},
    .window_size{3},
    .expected_labels{"CCC"s, "AGGG"s, "C"s, "GGG"s}
}));
