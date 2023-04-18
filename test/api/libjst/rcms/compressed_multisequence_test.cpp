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

#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include <libcontrib/seqan/alphabet.hpp>

#include <libjst/coverage/concept.hpp>
#include <libjst/coverage/bit_coverage.hpp>
#include <libjst/rcms/compressed_multisequence.hpp>

using jst::contrib::operator""_dna4;

struct compressed_multisequence_test : public ::testing::Test {
    using source_type = std::vector<jst::contrib::dna4>;
    using coverage_type = libjst::bit_coverage<uint32_t>;
    using coverage_domain_type = libjst::coverage_domain_t<coverage_type>;

    using test_type = libjst::compressed_multisequence<source_type, coverage_type>;
};

TEST_F(compressed_multisequence_test, range_concept) {
    EXPECT_TRUE(std::ranges::random_access_range<test_type>);
}

TEST_F(compressed_multisequence_test, construct) {
    source_type src{"AAAAAAAAAAAAAAA"_dna4};
    coverage_domain_type domain{0, 10};
    test_type multisequence{src, domain};

    EXPECT_EQ(multisequence.source(), source_type{"AAAAAAAAAAAAAAA"_dna4});
    EXPECT_TRUE(multisequence.coverage_domain() == domain);
}

TEST_F(compressed_multisequence_test, insert_snv) {
    source_type src{"AAAAAAAAAAAAAAA"_dna4};

    test_type multisequence{src, coverage_domain_type{0, 10}};

    using value_type = std::ranges::range_value_t<test_type>;
    coverage_type test_coverage{{0, 1, 2}, multisequence.coverage_domain()};

    auto it = multisequence.insert(value_type{libjst::breakpoint{3, 1}, "T"_dna4, test_coverage});
    EXPECT_EQ(libjst::low_breakend(*it), 3);
    EXPECT_EQ(libjst::high_breakend(*it), 4);
    EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), "T"_dna4));
    EXPECT_EQ(libjst::coverage(*it), test_coverage);
}

TEST_F(compressed_multisequence_test, insert_insertion) {
    source_type src{"AAAAAAAAAAAAAAA"_dna4};

    test_type multisequence{src, coverage_domain_type{0, 10}};

    using value_type = std::ranges::range_value_t<test_type>;
    coverage_type test_coverage{{0, 1, 2}, multisequence.coverage_domain()};

    auto it = multisequence.insert(value_type{libjst::breakpoint{3, 0}, "TCGT"_dna4, test_coverage});
    EXPECT_EQ(libjst::low_breakend(*it), 3);
    EXPECT_EQ(libjst::high_breakend(*it), 3);
    EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), "TCGT"_dna4));
    EXPECT_EQ(libjst::coverage(*it), test_coverage);
}

TEST_F(compressed_multisequence_test, insert_deletion) {
    source_type src{"AAAAAAAAAAAAAAA"_dna4};

    test_type multisequence{src, coverage_domain_type{0, 10}};

    using value_type = std::ranges::range_value_t<test_type>;
    coverage_type test_coverage{{0, 1, 2}, multisequence.coverage_domain()};

    auto it = multisequence.insert(value_type{libjst::breakpoint{3, 3}, ""_dna4, test_coverage});
    EXPECT_EQ(libjst::low_breakend(*it), 3);
    EXPECT_EQ(libjst::high_breakend(*it), 6);
    EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), ""_dna4));
    EXPECT_EQ(libjst::coverage(*it), test_coverage);
}

TEST_F(compressed_multisequence_test, iterate) {
    source_type src{"AAAAAAAAAAAAAAA"_dna4};
    using value_type = std::ranges::range_value_t<test_type>;
    coverage_domain_type domain{0, 10};
    coverage_type full_coverage{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, domain};

    { // empty
        test_type rcms{src, domain};
        auto it = rcms.begin();
        EXPECT_EQ(libjst::low_breakend(*it), 0);
        EXPECT_EQ(libjst::high_breakend(*it), 0);
        EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), ""_dna4));
        EXPECT_EQ(libjst::coverage(*it), full_coverage);

        ++it;
        EXPECT_EQ(libjst::low_breakend(*it), std::ranges::size(src));
        EXPECT_EQ(libjst::high_breakend(*it), std::ranges::size(src));
        EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), ""_dna4));
        EXPECT_EQ(libjst::coverage(*it), full_coverage);
        ++it;
        EXPECT_TRUE(it == rcms.end());
    }

    { // multiple variants
        test_type rcms{src, domain};
        coverage_type test_coverage{{0, 1, 2}, rcms.coverage_domain()};
        rcms.insert(value_type{libjst::breakpoint{9, 1}, "T"_dna4, test_coverage});
        rcms.insert(value_type{libjst::breakpoint{5, 1}, "T"_dna4, test_coverage});
        rcms.insert(value_type{libjst::breakpoint{1, 1}, "T"_dna4, test_coverage});
        rcms.insert(value_type{libjst::breakpoint{3, 1}, "T"_dna4, test_coverage});

        auto it = rcms.begin();
        EXPECT_EQ(libjst::low_breakend(*it), 0);
        EXPECT_EQ(libjst::high_breakend(*it), 0);
        EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), ""_dna4));
        EXPECT_EQ(libjst::coverage(*it), full_coverage);

        ++it;
        EXPECT_EQ(libjst::low_breakend(*it), 1);
        EXPECT_EQ(libjst::high_breakend(*it), 2);
        EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), "T"_dna4));
        EXPECT_EQ(libjst::coverage(*it), test_coverage);

        ++it;
        EXPECT_EQ(libjst::low_breakend(*it), 3);
        EXPECT_EQ(libjst::high_breakend(*it), 4);
        EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), "T"_dna4));
        EXPECT_EQ(libjst::coverage(*it), test_coverage);

        ++it;
        EXPECT_EQ(libjst::low_breakend(*it), 5);
        EXPECT_EQ(libjst::high_breakend(*it), 6);
        EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), "T"_dna4));
        EXPECT_EQ(libjst::coverage(*it), test_coverage);

        ++it;
        EXPECT_EQ(libjst::low_breakend(*it), 9);
        EXPECT_EQ(libjst::high_breakend(*it), 10);
        EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), "T"_dna4));
        EXPECT_EQ(libjst::coverage(*it), test_coverage);

        ++it;
        EXPECT_EQ(libjst::low_breakend(*it), std::ranges::size(src));
        EXPECT_EQ(libjst::high_breakend(*it), std::ranges::size(src));
        EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), ""_dna4));
        EXPECT_EQ(libjst::coverage(*it), full_coverage);
        ++it;
        EXPECT_TRUE(it == rcms.end());
    }
}

TEST_F(compressed_multisequence_test, source) {
    source_type src{"AAAAAAAAAAAAAAA"_dna4};

    test_type rcms{};
    EXPECT_RANGE_EQ(rcms.source(), ""_dna4);

    rcms = test_type{src, coverage_domain_type{0, 10}};
    EXPECT_RANGE_EQ(rcms.source(), src);
}

TEST_F(compressed_multisequence_test, empty) {
    source_type src{"AAAAAAAAAAAAAAA"_dna4};
    using value_type = std::ranges::range_value_t<test_type>;

    test_type rcms{};
    EXPECT_TRUE(std::ranges::empty(rcms));

    rcms = test_type{src, coverage_domain_type{0, 10}};
    EXPECT_FALSE(std::ranges::empty(rcms));

    coverage_type test_coverage{{0, 1, 2}, rcms.coverage_domain()};
    rcms.insert(value_type{libjst::breakpoint{9, 1}, "T"_dna4, test_coverage});

    EXPECT_FALSE(std::ranges::empty(rcms));
}

TEST_F(compressed_multisequence_test, size) {
    source_type src{"AAAAAAAAAAAAAAA"_dna4};
    using value_type = std::ranges::range_value_t<test_type>;

    test_type rcms{};
    EXPECT_EQ(rcms.size(), 0u);

    rcms = test_type{src, coverage_domain_type{0, 10}};
    EXPECT_EQ(rcms.size(), 2u);

    coverage_type test_coverage{{0, 1, 2}, rcms.coverage_domain()};
    rcms.insert(value_type{libjst::breakpoint{9, 1}, "T"_dna4, test_coverage});
    EXPECT_EQ(rcms.size(), 3u);
    rcms.insert(value_type{libjst::breakpoint{5, 1}, "T"_dna4, test_coverage});
    rcms.insert(value_type{libjst::breakpoint{1, 1}, "T"_dna4, test_coverage});
    rcms.insert(value_type{libjst::breakpoint{3, 1}, "T"_dna4, test_coverage});
    EXPECT_EQ(rcms.size(), 6u);
}

TEST_F(compressed_multisequence_test, has_conflicts) {
    source_type src{"AAAAAAAAAAAAAAA"_dna4};
    using value_type = std::ranges::range_value_t<test_type>;

    test_type rcms{src, coverage_domain_type{0, 10}};
    coverage_type test_coverage{{0, 1, 2}, rcms.coverage_domain()};
    rcms.insert(value_type{libjst::breakpoint{9, 1}, "T"_dna4, test_coverage});
    rcms.insert(value_type{libjst::breakpoint{5, 1}, "T"_dna4, test_coverage});
    rcms.insert(value_type{libjst::breakpoint{1, 1}, "T"_dna4, test_coverage});
    rcms.insert(value_type{libjst::breakpoint{3, 1}, "T"_dna4, test_coverage});

    EXPECT_FALSE((rcms.has_conflicts(value_type{libjst::breakpoint{0, 1}, "T"_dna4, test_coverage})));
    EXPECT_FALSE((rcms.has_conflicts(value_type{libjst::breakpoint{14, 1}, "T"_dna4, test_coverage})));
    EXPECT_FALSE((rcms.has_conflicts(value_type{libjst::breakpoint{2, 1}, "T"_dna4, test_coverage})));
    EXPECT_FALSE((rcms.has_conflicts(value_type{libjst::breakpoint{1, 1}, "T"_dna4, coverage_type{{3, 4, 5}, rcms.coverage_domain()}})));
    EXPECT_TRUE((rcms.has_conflicts(value_type{libjst::breakpoint{1, 1}, "T"_dna4, coverage_type{{1, 3, 4, 5}, rcms.coverage_domain()}})));
    EXPECT_TRUE((rcms.has_conflicts(value_type{libjst::breakpoint{9, 1}, "T"_dna4, coverage_type{{2, 9}, rcms.coverage_domain()}})));
}

TEST_F(compressed_multisequence_test, serialise) {
    using value_type = std::ranges::range_value_t<test_type>;

    source_type src{"AAAAAAAAAAAAAAA"_dna4};
    coverage_domain_type domain{0, 10};
    coverage_type full_coverage{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, domain};

    test_type rcms_out{src, coverage_domain_type{0, 10}};
    coverage_type test_coverage{{0, 1, 2}, rcms_out.coverage_domain()};
    rcms_out.insert(value_type{libjst::breakpoint{9, 1}, "T"_dna4, test_coverage});
    rcms_out.insert(value_type{libjst::breakpoint{5, 1}, "T"_dna4, test_coverage});
    rcms_out.insert(value_type{libjst::breakpoint{1, 1}, "T"_dna4, test_coverage});
    rcms_out.insert(value_type{libjst::breakpoint{3, 1}, "T"_dna4, test_coverage});

    std::stringstream buffer{};
    { // writing to output string stream
        cereal::JSONOutputArchive oarch{buffer};
        rcms_out.save(oarch);
    }

    test_type rcms_in{};
    { // writing to output string stream
        cereal::JSONInputArchive iarch{buffer};
        rcms_in.load(iarch);
    }

    auto it = rcms_in.begin();

    EXPECT_EQ(libjst::low_breakend(*it), 0);
    EXPECT_EQ(libjst::high_breakend(*it), 0);
    EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), ""_dna4));
    EXPECT_EQ(libjst::coverage(*it), full_coverage);

    ++it;
    EXPECT_EQ(libjst::low_breakend(*it), 1);
    EXPECT_EQ(libjst::high_breakend(*it), 2);
    EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), "T"_dna4));
    EXPECT_EQ(libjst::coverage(*it), test_coverage);

    ++it;
    EXPECT_EQ(libjst::low_breakend(*it), 3);
    EXPECT_EQ(libjst::high_breakend(*it), 4);
    EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), "T"_dna4));
    EXPECT_EQ(libjst::coverage(*it), test_coverage);

    ++it;
    EXPECT_EQ(libjst::low_breakend(*it), 5);
    EXPECT_EQ(libjst::high_breakend(*it), 6);
    EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), "T"_dna4));
    EXPECT_EQ(libjst::coverage(*it), test_coverage);

    ++it;
    EXPECT_EQ(libjst::low_breakend(*it), 9);
    EXPECT_EQ(libjst::high_breakend(*it), 10);
    EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), "T"_dna4));
    EXPECT_EQ(libjst::coverage(*it), test_coverage);

    ++it;
    EXPECT_EQ(libjst::low_breakend(*it), std::ranges::size(src));
    EXPECT_EQ(libjst::high_breakend(*it), std::ranges::size(src));
    EXPECT_TRUE(std::ranges::equal(libjst::alt_sequence(*it), ""_dna4));
    EXPECT_EQ(libjst::coverage(*it), full_coverage);

    ++it;
    EXPECT_TRUE(it == rcms_in.end());
}
