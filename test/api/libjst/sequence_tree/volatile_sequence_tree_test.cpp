// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <string>

#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>

#include "rcs_store_mock.hpp"

namespace jst::test::volatile_sequence_tree {

using source_t = std::string;
using variant_t = jst::test::variant<int, source_t, int, libjst::bit_vector<>>;

struct expected_root {
    source_t sequence{};
};

struct fixture {
    source_t source{};
    std::vector<variant_t> variants{};
    std::size_t coverage_size{};
    std::size_t window_size{};

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

} // namespace jst::test::volatile_sequence_tree

using namespace std::literals;

using fixture = jst::test::volatile_sequence_tree::fixture;
using variant_t = jst::test::volatile_sequence_tree::variant_t;
struct volatile_sequence_tree_test : public jst::test::volatile_sequence_tree::test
{
    using typename jst::test::volatile_sequence_tree::test::rcs_store_t;
    using jst::test::volatile_sequence_tree::test::get_mock;
    using jst::test::volatile_sequence_tree::test::GetParam;

    auto make_tree() const noexcept {
        auto const & rcs_mock = get_mock();
        auto mock_tree = libjst::volatile_tree{rcs_mock} | libjst::labelled() | libjst::coloured();
        return mock_tree;
    }

    auto expected_root_label() const noexcept {
        // find first variant
        int32_t label_end = std::ranges::ssize(get_mock().source());
        if (!GetParam().variants.empty()) {
            label_end = libjst::left_breakpoint(*GetParam().variants.begin());
        }

        return std::ranges::subrange{_mock.source().begin(),
                                     std::ranges::next(_mock.source().begin(), label_end, _mock.source().end())};
    }

    auto expected_reference_path() const noexcept {
        return _mock.source();
    }
};

// ----------------------------------------------------------------------------
// Test case definitions
// ----------------------------------------------------------------------------

TEST_P(volatile_sequence_tree_test, root_label) {
    auto tree = make_tree();
    auto r = tree.root();
    auto s = tree.sink();
    EXPECT_TRUE(r != s);

    EXPECT_TRUE(std::ranges::equal(r.label(), expected_root_label()));
}

TEST_P(volatile_sequence_tree_test, reference_path_label) {
    auto tree = make_tree();
    auto v = tree.root();
    auto s = tree.sink();
    EXPECT_TRUE(v != s);

    std::string actual_ref_path{};

    auto append_label = [&] (auto const & w) {
        auto lbl = w.label();
        actual_ref_path.append(std::string{lbl.begin(), lbl.end()});
    };

    while (v != s) {
        append_label(v);
        v = *v.next_ref();
    }

    EXPECT_TRUE(v == s);
    EXPECT_TRUE(std::ranges::equal(actual_ref_path, expected_reference_path()));
}

TEST_P(volatile_sequence_tree_test, reference_path_coverage) {
    using coverage_t = libjst::variant_coverage_t<variant_t>;
    auto tree = make_tree();
    auto v = tree.root();
    auto s = tree.sink();
    EXPECT_TRUE(v != s);

    coverage_t base_cov(GetParam().coverage_size, true);
    while (v != s) {
        EXPECT_TRUE(std::ranges::equal(v.coverage(), base_cov));
        v = *v.next_ref();
    }

    EXPECT_TRUE(v == s);
}

// ----------------------------------------------------------------------------
// Test values
// ----------------------------------------------------------------------------

INSTANTIATE_TEST_SUITE_P(no_variant, volatile_sequence_tree_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{},
    .coverage_size{4},
    .window_size{4}
}));

INSTANTIATE_TEST_SUITE_P(snv_first_base, volatile_sequence_tree_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{
        variant_t{.position{0}, .insertion{"x"s}, .deletion{1}, .coverage{1, 0, 0, 0}}
    },
    .coverage_size{4},
    .window_size{4}
}));

INSTANTIATE_TEST_SUITE_P(snv_last_base, volatile_sequence_tree_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{
        variant_t{.position{7}, .insertion{"x"s}, .deletion{1}, .coverage{1, 0, 0, 0}}
    },
    .coverage_size{4},
    .window_size{4}
}));

INSTANTIATE_TEST_SUITE_P(snv_middle, volatile_sequence_tree_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{
        variant_t{.position{4}, .insertion{"x"s}, .deletion{1}, .coverage{1, 0, 0, 0}}
    },
    .coverage_size{4},
    .window_size{4}
}));
