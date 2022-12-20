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

struct naive_matcher {
    source_t needle;

    template <typename seq_t, typename callback_t>
    constexpr void operator()(seq_t && seq, callback_t && callback) const {
        seqan3::debug_stream << "Hystk: " << seq << "\n";
        auto it = std::ranges::next(std::ranges::begin(seq), window_size(), std::ranges::end(seq));
        for (; it != std::ranges::end(seq); ++it) {
            if (auto [has_hit, pos] = find_impl(it); has_hit)
                callback(pos);
        }
    }

    std::size_t window_size() const noexcept {
        return needle.size();
    }

private:

    template <typename it_t>
    constexpr auto find_impl(it_t it) const noexcept {
        auto hst_it = std::make_reverse_iterator(std::move(it));
        auto rev_needle = needle | std::views::reverse;
        for (auto ndl_it = std::ranges::begin(rev_needle); ndl_it != std::ranges::end(rev_needle); ++ndl_it, ++hst_it) {
            if (*ndl_it != *hst_it) return std::pair{false, hst_it.base()};
        }
        return std::pair{true, hst_it.base()};
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
        return jst::test::polymorphic_sequence_searcher::naive_matcher{GetParam().needle};
    }
};

// ----------------------------------------------------------------------------
// Test case definitions
// ----------------------------------------------------------------------------

TEST_P(polymorphic_sequence_searcher_test, search) {
    auto searcher = make_searcher();
    auto pattern = get_pattern();

    std::vector<std::size_t> actual_occurrences{};
    searcher(pattern, [&] (auto && lbl_it, auto const & cargo) {
        actual_occurrences.push_back((lbl_it - cargo.sequence().begin()));
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
    .expected_occurrences{2}
}));

INSTANTIATE_TEST_SUITE_P(single_snv_variant, polymorphic_sequence_searcher_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{variant_t{.position{4}, .insertion{"O"}, .deletion{1}, .coverage{1,1,0,0}}},
    .coverage_size{4},
    .needle{"aaOb"},
    .expected_occurrences{1}
}));
