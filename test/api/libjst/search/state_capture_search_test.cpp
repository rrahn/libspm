// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <concepts>
#include <string>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <libjst/matcher/shiftor_matcher_restorable.hpp>
#include <libjst/search/polymorphic_sequence_searcher.hpp>

#include "../mock/rcs_store_mock.hpp"

namespace jst::test::polymorphic_sequence_searcher {

using source_t = std::string;
using variant_t = jst::test::variant<libjst::breakpoint, source_t, int, libjst::bit_vector<>>;

struct fixture {
    source_t source{};
    std::vector<variant_t> variants{};
    std::size_t coverage_size{};
    source_t needle{};
    std::vector<std::size_t> expected_occurrences{};

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

    virtual void SetUp() override {
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

} // namespace jst::test::polymorphic_sequence_searcher

using namespace std::literals;

using fixture = jst::test::polymorphic_sequence_searcher::fixture;
using variant_t = jst::test::polymorphic_sequence_searcher::variant_t;

struct polymorphic_sequence_searcher_test : public jst::test::polymorphic_sequence_searcher::test
{
    using base_test_t = jst::test::polymorphic_sequence_searcher::test;

    using typename base_test_t::rcs_store_t;
    using base_test_t::get_mock;
    using base_test_t::GetParam;

    virtual void SetUp() override {
        base_test_t::SetUp();
    }

    auto make_searcher() const noexcept {
        auto const & rcs_mock = get_mock();
        return libjst::polymorphic_sequence_searcher{rcs_mock};
    }

    auto get_pattern() const noexcept {
        return libjst::restorable_shiftor_matcher{GetParam().needle};
    }
};

// ----------------------------------------------------------------------------
// Test case definitions
// ----------------------------------------------------------------------------

TEST_P(polymorphic_sequence_searcher_test, search) {
    auto searcher = make_searcher();
    auto pattern = get_pattern();

    EXPECT_TRUE(libjst::restorable_matcher<std::remove_reference_t<decltype(pattern)> &>);

    std::vector<std::size_t> actual_occurrences{};
    searcher(pattern, [&] (auto && lbl_it, [[maybe_unused]] auto const & cargo) {
        actual_occurrences.push_back(seqan::endPosition(lbl_it));
    });

    std::ranges::sort(actual_occurrences);
    EXPECT_RANGE_EQ(actual_occurrences, GetParam().expected_occurrences);
}

// ----------------------------------------------------------------------------
// Test values
// ----------------------------------------------------------------------------

INSTANTIATE_TEST_SUITE_P(no_variant, polymorphic_sequence_searcher_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{},
    .coverage_size{4},
    .needle{"aabb"},
    .expected_occurrences{6}
}));

INSTANTIATE_TEST_SUITE_P(single_snv_variant, polymorphic_sequence_searcher_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{variant_t{.position{4}, .insertion{"O"}, .deletion{1}, .coverage{1,1,0,0}}},
    .coverage_size{4},
    .needle{"aaOb"},
    .expected_occurrences{2}
}));

INSTANTIATE_TEST_SUITE_P(single_snv_variant_at_begin, polymorphic_sequence_searcher_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{variant_t{.position{0}, .insertion{"O"}, .deletion{1}, .coverage{1,1,0,0}}},
    .coverage_size{4},
    .needle{"Oaaa"},
    .expected_occurrences{4}
}));

INSTANTIATE_TEST_SUITE_P(single_snv_variant_at_end, polymorphic_sequence_searcher_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{variant_t{.position{7}, .insertion{"O"}, .deletion{1}, .coverage{1,1,0,0}}},
    .coverage_size{4},
    .needle{"bbbO"},
    .expected_occurrences{1}
}));

INSTANTIATE_TEST_SUITE_P(two_snv_variants_on_different_subtrees, polymorphic_sequence_searcher_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{variant_t{.position{1}, .insertion{"I"}, .deletion{1}, .coverage{1,1,0,0}},
              variant_t{.position{5}, .insertion{"J"}, .deletion{1}, .coverage{1,1,0,0}}},
    .coverage_size{4},
    .needle{"Iaab"},
    .expected_occurrences{4}
}));

INSTANTIATE_TEST_SUITE_P(two_snv_variants_on_same_subtree, polymorphic_sequence_searcher_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{variant_t{.position{1}, .insertion{"I"}, .deletion{1}, .coverage{1,1,0,0}},
              variant_t{.position{4}, .insertion{"J"}, .deletion{1}, .coverage{1,0,0,0}}},
    .coverage_size{4},
    .needle{"IaaJ"},
    .expected_occurrences{1}
}));

INSTANTIATE_TEST_SUITE_P(two_snv_variants_behind_each_other, polymorphic_sequence_searcher_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{variant_t{.position{3}, .insertion{"I"}, .deletion{1}, .coverage{1,1,0,0}},
              variant_t{.position{4}, .insertion{"J"}, .deletion{1}, .coverage{1,0,0,0}}},
    .coverage_size{4},
    .needle{"aIJb"},
    .expected_occurrences{2}
}));

INSTANTIATE_TEST_SUITE_P(two_snv_variants_mutual_exclusive, polymorphic_sequence_searcher_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{variant_t{.position{3}, .insertion{"I"}, .deletion{1}, .coverage{1,1,0,0}},
              variant_t{.position{4}, .insertion{"J"}, .deletion{1}, .coverage{0,0,1,1}}},
    .coverage_size{4},
    .needle{"aIbb"},
    .expected_occurrences{2}
}));

INSTANTIATE_TEST_SUITE_P(two_snv_variants_mutual_exclusive_at_same_position, polymorphic_sequence_searcher_test,
testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{variant_t{.position{4}, .insertion{"I"}, .deletion{1}, .coverage{1,1,0,0}},
              variant_t{.position{4}, .insertion{"J"}, .deletion{1}, .coverage{0,0,1,1}}},
    .coverage_size{4},
    .needle{"Jbbb"},
    .expected_occurrences{4}
}));

INSTANTIATE_TEST_SUITE_P(three_snv_variants_in_same_subtree, polymorphic_sequence_searcher_test,
testing::Values(fixture{
        //   01234567
    .source{"aaaabbbb"},
    .variants{variant_t{.position{3}, .insertion{"I"}, .deletion{1}, .coverage{1,1,0,0}},
              variant_t{.position{4}, .insertion{"J"}, .deletion{1}, .coverage{0,1,1,0}},
              variant_t{.position{5}, .insertion{"K"}, .deletion{1}, .coverage{0,1,0,1}}},
    .coverage_size{4},
    .needle{"aIJKb"},
    .expected_occurrences{2}
}));

INSTANTIATE_TEST_SUITE_P(three_snv_variants_in_same_subtree_two_on_same_position, polymorphic_sequence_searcher_test,
testing::Values(fixture{
        //   01234567
    .source{"aaaabbbb"},
    .variants{variant_t{.position{3}, .insertion{"I"}, .deletion{1}, .coverage{1,1,0,0}},
              variant_t{.position{5}, .insertion{"J"}, .deletion{1}, .coverage{1,0,1,0}},
              variant_t{.position{5}, .insertion{"K"}, .deletion{1}, .coverage{0,1,0,1}}},
    .coverage_size{4},
    .needle{"aIbKb"},
    .expected_occurrences{2}
}));
