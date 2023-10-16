// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>
#include <sstream>

#include <libjst/coverage/concept.hpp>
#include <libjst/coverage/bit_coverage.hpp>
#include <libjst/rcms/dna_compressed_multisequence.hpp>
#include <libjst/rcms/compressed_multisequence_reversed.hpp>

using namespace std::literals;

struct compressed_multisequence_reversed_test : public ::testing::Test {
    using source_type = std::string;
    using coverage_type = libjst::bit_coverage<uint32_t>;
    using coverage_domain_type = libjst::coverage_domain_t<coverage_type>;

    using wrapped_test_type = libjst::dna_compressed_multisequence<source_type, coverage_type>;
    using test_type = libjst::compressed_multisequence_reversed<wrapped_test_type>;
};

TEST_F(compressed_multisequence_reversed_test, range_concept) {
    EXPECT_TRUE(std::ranges::random_access_range<test_type>);
}

TEST_F(compressed_multisequence_reversed_test, construct) {
    source_type src{"AAAAAAAAGGGGGGG"s};
    coverage_domain_type domain{0, 10};
    wrapped_test_type multisequence{src, domain};
    test_type reverse_rcms{multisequence};
    EXPECT_TRUE(std::ranges::equal(reverse_rcms.source(), source_type{"GGGGGGGAAAAAAAA"s}));
    EXPECT_TRUE(reverse_rcms.coverage_domain() == domain);
}

TEST_F(compressed_multisequence_reversed_test, iterate) {
    source_type src{"AAAAAAAAAAAAAAA"s};
    using value_type = std::ranges::range_value_t<wrapped_test_type>;
    coverage_domain_type domain{0, 10};
    coverage_type full_coverage{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, domain};

    { // empty
        wrapped_test_type rcms{src, domain};
        test_type reverse_rcms{rcms};

        auto it = reverse_rcms.begin();
        EXPECT_EQ(libjst::low_breakend(*it), 0);
        EXPECT_EQ(libjst::high_breakend(*it), 0);
        EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), ""s));
        EXPECT_EQ(libjst::coverage(*it), full_coverage);

        ++it;
        EXPECT_EQ(libjst::low_breakend(*it), std::ranges::size(src));
        EXPECT_EQ(libjst::high_breakend(*it), std::ranges::size(src));
        EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), ""s));
        EXPECT_EQ(libjst::coverage(*it), full_coverage);
        ++it;
        EXPECT_TRUE(it == reverse_rcms.end());
    }

    { // multiple variants
        wrapped_test_type rcms{src, domain};
        coverage_type test_coverage{{0, 1, 2}, rcms.coverage_domain()};
        rcms.insert(value_type{libjst::breakpoint{9, 1}, "T"s, test_coverage});
        rcms.insert(value_type{libjst::breakpoint{5, 1}, "C"s, test_coverage});
        rcms.insert(value_type{libjst::breakpoint{1, 1}, "G"s, test_coverage});
        rcms.insert(value_type{libjst::breakpoint{3, 1}, "T"s, test_coverage});

        test_type reverse_rcms{rcms};

        auto it = reverse_rcms.begin();
        EXPECT_EQ(libjst::low_breakend(*it), 0);
        EXPECT_EQ(libjst::high_breakend(*it), 0);
        EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), ""s));
        EXPECT_EQ(libjst::coverage(*it), full_coverage);

        ++it;
        EXPECT_EQ(libjst::low_breakend(*it), 5);
        EXPECT_EQ(libjst::high_breakend(*it), 6);
        EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), "T"s));
        EXPECT_EQ(libjst::coverage(*it), test_coverage);

        ++it;
        EXPECT_EQ(libjst::low_breakend(*it), 9);
        EXPECT_EQ(libjst::high_breakend(*it), 10);
        EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), "C"s));
        EXPECT_EQ(libjst::coverage(*it), test_coverage);

        ++it;
        EXPECT_EQ(libjst::low_breakend(*it), 11);
        EXPECT_EQ(libjst::high_breakend(*it), 12);
        EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), "T"s));
        EXPECT_EQ(libjst::coverage(*it), test_coverage);

        ++it;
        EXPECT_EQ(libjst::low_breakend(*it), 13);
        EXPECT_EQ(libjst::high_breakend(*it), 14);
        EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), "G"s));
        EXPECT_EQ(libjst::coverage(*it), test_coverage);

        ++it;
        EXPECT_EQ(libjst::low_breakend(*it), std::ranges::size(src));
        EXPECT_EQ(libjst::high_breakend(*it), std::ranges::size(src));
        EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), ""s));
        EXPECT_EQ(libjst::coverage(*it), full_coverage);
        ++it;
        EXPECT_TRUE(it == reverse_rcms.end());
    }
}

TEST_F(compressed_multisequence_reversed_test, source) {
    source_type src{"AACCGGTTAAACCCG"s};

    wrapped_test_type rcms{};
    test_type reverse_rcms{rcms};

    EXPECT_TRUE(std::ranges::equal(reverse_rcms.source(), ""s));

    rcms = wrapped_test_type{src, coverage_domain_type{0, 10}};
    EXPECT_TRUE(std::ranges::equal(reverse_rcms.source(), "GCCCAAATTGGCCAA"s));
}

TEST_F(compressed_multisequence_reversed_test, empty) {
    source_type src{"AACCGGTTAAACCCG"s};
    using original_value_type = std::ranges::range_value_t<wrapped_test_type>;

    wrapped_test_type rcms{};
    test_type reverse_rcms{rcms};
    EXPECT_TRUE(std::ranges::empty(reverse_rcms));

    rcms = wrapped_test_type{src, coverage_domain_type{0, 10}};
    EXPECT_FALSE(std::ranges::empty(reverse_rcms));

    coverage_type test_coverage{{0, 1, 2}, rcms.coverage_domain()};
    rcms.insert(original_value_type{libjst::breakpoint{9, 1}, "T"s, test_coverage});

    EXPECT_FALSE(std::ranges::empty(reverse_rcms));
}

TEST_F(compressed_multisequence_reversed_test, size) {
    source_type src{"AACCGGTTAAACCCG"s};
    using original_value_type = std::ranges::range_value_t<wrapped_test_type>;

    wrapped_test_type rcms{};
    test_type reverse_rcms{rcms};
    EXPECT_EQ(reverse_rcms.size(), 0u);

    rcms = wrapped_test_type{src, coverage_domain_type{0, 10}};
    EXPECT_EQ(reverse_rcms.size(), 2u);

    coverage_type test_coverage{{0, 1, 2}, rcms.coverage_domain()};
    rcms.insert(original_value_type{libjst::breakpoint{9, 1}, "T"s, test_coverage});
    EXPECT_EQ(reverse_rcms.size(), 3u);
    rcms.insert(original_value_type{libjst::breakpoint{5, 1}, "T"s, test_coverage});
    rcms.insert(original_value_type{libjst::breakpoint{1, 1}, "T"s, test_coverage});
    rcms.insert(original_value_type{libjst::breakpoint{3, 1}, "T"s, test_coverage});
    EXPECT_EQ(reverse_rcms.size(), 6u);
}
