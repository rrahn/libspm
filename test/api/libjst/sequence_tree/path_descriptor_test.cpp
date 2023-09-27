// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <libjst/sequence_tree/path_descriptor.hpp>

struct path_descriptor_test : public ::testing::Test {
};

TEST(path_descriptor_test, use_case) {
    libjst::alternate_path_descriptor descr{};

    EXPECT_EQ(descr.size(), 1);
    EXPECT_EQ(descr.max_size(), 256);
    for (unsigned i = 1; i < descr.max_size(); ++i) {
        descr.next();
        bool is_ref = (i & 1);
        if (is_ref) {
            descr.set_ref();
        } else {
            descr.set_alt();
        }
        EXPECT_EQ(descr.size(), i + 1);
    }
    EXPECT_EQ(descr.size(), descr.max_size());

    std::vector<bool> expected_path{};
    expected_path.resize(descr.max_size());
    bool v{true};
    for (auto && elem : expected_path) {
        elem = v;
        v = !v;
    }

    EXPECT_TRUE(std::ranges::equal(descr, expected_path));
}

// TEST(path_descriptor_test, )
