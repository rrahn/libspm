// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/../../../unit/range/iterator_test_template.hpp>

#include <libjst/journal.hpp>

// using segment_t = std::span<char>;
using key_type = uint32_t;
using mapped_type = libjst::detail::subrange_t<std::string &>;
using journal_type = libjst::journal<key_type, std::string &>;
using journal_sequence_type = decltype(std::declval<journal_type &>().sequence());
using journal_sequence_iterator = std::ranges::iterator_t<journal_sequence_type>;

template <>
struct iterator_fixture<journal_sequence_iterator> : public ::testing::Test
{
    std::string reference{"aaaaaaaa"};
    std::string ins_segment{"ccccgggggggg"};
    std::string repl_segment{"tttt"};

    using iterator_tag = std::random_access_iterator_tag;
    static constexpr bool const_iterable = false;

    journal_type journal{};
    journal_sequence_type test_range{};
    std::string expected_range{"aaaaccccggggtttt"};

    void SetUp() override
    {
        journal = journal_type{reference}; // aaaa
        journal.record_insertion(4, ins_segment); // aaaaccccggggggggaaaa
        journal.record_substitution(16, repl_segment); // aaaaccccggggggggtttt
        journal.record_deletion(9, 4); // aaaaccccggggtttt
        test_range = journal.sequence();

        EXPECT_RANGE_EQ(test_range, expected_range);
    }
};

INSTANTIATE_TYPED_TEST_SUITE_P(journal_sequence_iterator_test,
                               iterator_fixture,
                               ::testing::Types<journal_sequence_iterator>, );
