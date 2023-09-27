// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <libjst/rcms/packed_breakend_key.hpp>
#include <libjst/utility/multi_invocable.hpp>

struct packed_breakend_key : public ::testing::Test {
    using test_type = libjst::packed_breakend_key<>;
};

TEST_F(packed_breakend_key, construct_snv) {
    test_type test{static_cast<uint8_t>(2), 3000};

    EXPECT_FALSE(test.is_indel());
    EXPECT_EQ(test.snv_value(), 2);
    EXPECT_EQ(test.position(), 3000);
}

TEST_F(packed_breakend_key, construct_deletion_low) {
    test_type test{libjst::indel_breakend_kind::deletion_low, 1236};

    EXPECT_TRUE(test.is_indel());
    EXPECT_EQ(test.indel_kind(), libjst::indel_breakend_kind::deletion_low);
    EXPECT_EQ(test.position(), 1236);
}

TEST_F(packed_breakend_key, construct_deletion_high) {
    test_type test{libjst::indel_breakend_kind::deletion_high, 6321};

    EXPECT_TRUE(test.is_indel());
    EXPECT_EQ(test.indel_kind(), libjst::indel_breakend_kind::deletion_high);
    EXPECT_EQ(test.position(), 6321);
}

TEST_F(packed_breakend_key, construct_insertion) {
    test_type test{libjst::indel_breakend_kind::insertion_low, 0};

    EXPECT_TRUE(test.is_indel());
    EXPECT_EQ(test.indel_kind(), libjst::indel_breakend_kind::insertion_low);
    EXPECT_EQ(test.position(), 0);
}

TEST_F(packed_breakend_key, max_position) {
    using value_t = typename test_type::underlying_type;
    value_t max_position = (1 << 29) - 1;
    test_type test{static_cast<uint8_t>(0), max_position};
    EXPECT_EQ(test.position(), max_position);

    test_type test2{static_cast<uint8_t>(0), max_position + 1};
    EXPECT_EQ(test2.position(), 0);
}

TEST_F(packed_breakend_key, visit) {
    using value_t = typename test_type::underlying_type;
    test_type test_snv{static_cast<uint8_t>(2), 3000};
    test_snv.visit(libjst::multi_invocable{
        [] (libjst::indel_breakend_kind) { FAIL() << "Expected snv."; },
        [] (value_t snv) { EXPECT_EQ(snv, 2); }
    });

    test_type test_del{libjst::indel_breakend_kind::deletion_high, 6321};
    test_del.visit(libjst::multi_invocable{
        [] (libjst::indel_breakend_kind indel) { EXPECT_EQ(indel, libjst::indel_breakend_kind::deletion_high); },
        [] (value_t ) { FAIL() << "Expected indel."; }
    });
}

TEST_F(packed_breakend_key, equal) {
    test_type test_snv{static_cast<uint8_t>(2), 3000};
    test_type test_del{libjst::indel_breakend_kind::deletion_high, 6321};

    EXPECT_EQ(test_snv, test_snv);
    EXPECT_EQ(test_del, test_del);
}

TEST_F(packed_breakend_key, unequal) {
    test_type test_snv{static_cast<uint8_t>(2), 3000};
    test_type test_del{libjst::indel_breakend_kind::deletion_high, 6321};
    test_type test_ins{libjst::indel_breakend_kind::insertion_low, 0};

    EXPECT_NE(test_snv, test_del);
    EXPECT_NE(test_snv, test_ins);
    EXPECT_NE(test_del, test_ins);
}


TEST_F(packed_breakend_key, less) {
    test_type{static_cast<uint8_t>(2), 3000};
    test_type test_del{libjst::indel_breakend_kind::deletion_high, 6321};
    test_type{libjst::indel_breakend_kind::insertion_low, 0};

    EXPECT_LT((test_type{static_cast<uint8_t>(2), 3000}), (test_type{static_cast<uint8_t>(3), 3000}));
    EXPECT_LT((test_type{static_cast<uint8_t>(2), 3000}), (test_type{libjst::indel_breakend_kind::insertion_low, 3001}));
    EXPECT_LT((test_type{static_cast<uint8_t>(2), 3000}), (test_type{libjst::indel_breakend_kind::deletion_low, 3000}));
    EXPECT_LT((test_type{static_cast<uint8_t>(2), 3000}), (test_type{libjst::indel_breakend_kind::deletion_high, 3001}));
    EXPECT_LT((test_type{libjst::indel_breakend_kind::insertion_low, 3000}),
              (test_type{static_cast<uint8_t>(0), 3000}));
    EXPECT_LT((test_type{libjst::indel_breakend_kind::insertion_low, 3000}),
              (test_type{libjst::indel_breakend_kind::insertion_low, 3001}));
    EXPECT_LT((test_type{libjst::indel_breakend_kind::insertion_low, 3000}),
              (test_type{libjst::indel_breakend_kind::deletion_low, 3000}));
    EXPECT_LT((test_type{libjst::indel_breakend_kind::insertion_low, 3000}),
              (test_type{libjst::indel_breakend_kind::deletion_high, 3001}));
    EXPECT_LT((test_type{libjst::indel_breakend_kind::deletion_low, 3000}),
              (test_type{libjst::indel_breakend_kind::deletion_low, 3001}));
    EXPECT_LT((test_type{libjst::indel_breakend_kind::deletion_low, 3000}),
              (test_type{libjst::indel_breakend_kind::deletion_high, 3001}));
    EXPECT_LT((test_type{libjst::indel_breakend_kind::deletion_high, 3000}),
              (test_type{libjst::indel_breakend_kind::deletion_high, 3001}));
    EXPECT_LT((test_type{libjst::indel_breakend_kind::deletion_high, 3000}),
              (test_type{libjst::indel_breakend_kind::deletion_low, 3000}));
    EXPECT_LT((test_type{libjst::indel_breakend_kind::deletion_high, 3000}),
              (test_type{libjst::indel_breakend_kind::insertion_low, 3000}));
    EXPECT_LT((test_type{libjst::indel_breakend_kind::deletion_high, 3000}),
              (test_type{static_cast<uint8_t>(0), 3000}));

}
