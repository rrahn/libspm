// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <cereal/archives/json.hpp> // output archive for testing

#include <seqan3/alphabet/adaptation/char.hpp>

#include <libjst/detail/delta_event_shared.hpp>

struct delta_event_shared_fixture : public ::testing::Test
{
    using delta_event_shared_t = libjst::detail::delta_event_shared<char>;
    using alphabet_t = typename delta_event_shared_t::alphabet_type;
    using insertion_t = typename delta_event_shared_t::insertion_type;
    using substitution_t = typename delta_event_shared_t::substitution_type;
    using deletion_t = typename delta_event_shared_t::deletion_type;
    using delta_event_t = typename delta_event_shared_t::delta_event_type;
    using coverage_t = typename delta_event_shared_t::coverage_type;
    using position_t = typename delta_event_shared_t::position_type;

    static constexpr std::string_view expected_substitution_archive =
R"json({
    "value0": {
        "value0": 23,
        "value1": {
            "index": 2,
            "data": {
                "value0": {
                    "value0": [
                        97,
                        98,
                        99,
                        100
                    ]
                }
            }
        }
    },
    "value1": {
        "value0": [
            5
        ],
        "value1": 4
    }
})json";

    static constexpr std::string_view expected_insertion_archive =
R"json({
    "value0": {
        "value0": 5,
        "value1": {
            "index": 0,
            "data": {
                "value0": {
                    "value0": [
                        105,
                        106,
                        107,
                        108,
                        109
                    ]
                }
            }
        }
    },
    "value1": {
        "value0": [
            5
        ],
        "value1": 4
    }
})json";

    static constexpr std::string_view expected_deletion_archive =
R"json({
    "value0": {
        "value0": 100,
        "value1": {
            "index": 3,
            "data": {
                "value0": {
                    "value0": 10
                }
            }
        }
    },
    "value1": {
        "value0": [
            5
        ],
        "value1": 4
    }
})json";
};

TEST_F(delta_event_shared_fixture, construction)
{
    EXPECT_TRUE(std::is_default_constructible_v<delta_event_shared_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<delta_event_shared_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<delta_event_shared_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<delta_event_shared_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<delta_event_shared_t>);
    EXPECT_TRUE(std::is_nothrow_destructible_v<delta_event_shared_t>);

    // Special constructor.
    EXPECT_TRUE((std::is_constructible_v<delta_event_shared_t, delta_event_t, coverage_t>));
    EXPECT_TRUE((std::is_constructible_v<delta_event_shared_t, position_t, insertion_t, coverage_t>));
    EXPECT_TRUE((std::is_constructible_v<delta_event_shared_t, position_t, deletion_t, coverage_t>));
    EXPECT_TRUE((std::is_constructible_v<delta_event_shared_t, position_t, substitution_t, coverage_t>));
}

TEST_F(delta_event_shared_fixture, construct_from_substitution)
{
    using namespace std::literals;

    delta_event_shared_t node1{delta_event_t{position_t{.offset = 10u},
                                             substitution_t{"abc"s}},
                                             coverage_t{0, 1, 1, 0, 0, 1, 1}};

    EXPECT_EQ(node1.position(), position_t{.offset = 10u});
    EXPECT_TRUE(node1.is_substitution());

    delta_event_shared_t node2{position_t{.offset = 10u}, substitution_t{"abc"s}, coverage_t{0, 1, 1, 0, 0, 1, 1}};
    EXPECT_EQ(node2.position(), position_t{.offset = 10u});
    EXPECT_TRUE(node2.is_substitution());
}

TEST_F(delta_event_shared_fixture, construct_from_insertion)
{
    using namespace std::literals;

    delta_event_shared_t node1{delta_event_t{position_t{.offset = 10u},
                               insertion_t{"abc"s}},
                               coverage_t{0, 1, 1, 0, 0, 1, 1}};

    EXPECT_EQ(node1.position(), position_t{.offset = 10u});
    EXPECT_TRUE(node1.is_insertion());

    delta_event_shared_t node2{position_t{.offset = 10u}, insertion_t{"abc"s}, coverage_t{0, 1, 1, 0, 0, 1, 1}};

    EXPECT_EQ(node2.position(), position_t{.offset = 10u});
    EXPECT_TRUE(node2.is_insertion());
}

TEST_F(delta_event_shared_fixture, construct_from_deletion)
{
    delta_event_shared_t node1{delta_event_t{position_t{.offset = 10u},
                                             deletion_t{1}},
                                             coverage_t{0, 1, 1, 0, 0, 1, 1}};

    EXPECT_EQ(node1.position(), position_t{.offset = 10u});
    EXPECT_TRUE(node1.is_deletion());

    delta_event_shared_t node2{position_t{.offset = 10u}, deletion_t{1}, coverage_t{0, 1, 1, 0, 0, 1, 1}};

    EXPECT_EQ(node2.position(), position_t{.offset = 10u});
    EXPECT_TRUE(node2.is_deletion());
}

TEST_F(delta_event_shared_fixture, position)
{
    delta_event_shared_t node{delta_event_t{position_t{.offset = 10u},
                                            deletion_t{1}},
                                            coverage_t{0, 1, 1, 0, 0, 1, 1}};

    EXPECT_EQ(node.position(), position_t{.offset = 10u});
    EXPECT_EQ(std::as_const(node).position(), position_t{.offset = 10u});
}

TEST_F(delta_event_shared_fixture, coverage)
{
    delta_event_shared_t node{delta_event_t{position_t{.offset = 10u},
                                            deletion_t{1}},
                                            coverage_t{0, 1, 1, 0, 0, 1, 1}};

    EXPECT_EQ(node.coverage(), (coverage_t{0, 1, 1, 0, 0, 1, 1}));
    EXPECT_EQ(std::as_const(node).coverage(), (coverage_t{0, 1, 1, 0, 0, 1, 1}));
}

TEST_F(delta_event_shared_fixture, equality)
{
    using namespace std::literals;

    delta_event_shared_t node1{delta_event_t{position_t{.offset = 10u},
                                             deletion_t{1}},
                                             coverage_t{0, 1, 1, 0, 0, 1, 1}};
    delta_event_shared_t node2{delta_event_t{position_t{.offset = 10u},
                                             deletion_t{1}},
                                             coverage_t{0, 1, 1, 0, 0, 1, 0}};
    delta_event_shared_t node3{delta_event_t{position_t{.offset = 9u},
                                             deletion_t{1}},
                                             coverage_t{0, 1, 1, 0, 0, 1, 0}};
    delta_event_shared_t node4{delta_event_t{position_t{.offset = 9u},
                                             substitution_t{"a"s}},
                                             coverage_t{0, 1, 1, 0, 0, 1, 0}};

    EXPECT_EQ(std::as_const(node1), std::as_const(node1));
    EXPECT_NE(std::as_const(node1), std::as_const(node2));
    EXPECT_NE(std::as_const(node1), std::as_const(node3));
    EXPECT_NE(std::as_const(node1), std::as_const(node4));

    EXPECT_NE(std::as_const(node2), std::as_const(node1));
    EXPECT_EQ(std::as_const(node2), std::as_const(node2));
    EXPECT_NE(std::as_const(node2), std::as_const(node3));
    EXPECT_NE(std::as_const(node2), std::as_const(node4));

    EXPECT_NE(std::as_const(node3), std::as_const(node1));
    EXPECT_NE(std::as_const(node3), std::as_const(node2));
    EXPECT_EQ(std::as_const(node3), std::as_const(node3));
    EXPECT_NE(std::as_const(node3), std::as_const(node4));

    EXPECT_NE(std::as_const(node4), std::as_const(node1));
    EXPECT_NE(std::as_const(node4), std::as_const(node2));
    EXPECT_NE(std::as_const(node4), std::as_const(node3));
    EXPECT_EQ(std::as_const(node4), std::as_const(node4));
}

TEST_F(delta_event_shared_fixture, less)
{
    using namespace std::literals;

    delta_event_shared_t node1{delta_event_t{position_t{.offset = 0u},
                                             deletion_t{2}},
                                             coverage_t{0, 1, 1, 0, 0, 1, 1}};
    delta_event_shared_t node2{delta_event_t{position_t{.offset = 0u},
                                             deletion_t{1}},
                                             coverage_t{0, 1, 1, 0, 0, 1, 1}};
    delta_event_shared_t node3{delta_event_t{position_t{.offset = 0u},
                                             deletion_t{1}},
                                             coverage_t{0, 1, 1, 0, 0, 1, 0}};
    delta_event_shared_t node4{delta_event_t{position_t{.offset = 0u},
                                             substitution_t{"abc"s}},
                                             coverage_t{0, 1, 1, 0, 0, 1, 1}};
    delta_event_shared_t node5{delta_event_t{position_t{.offset = 10u},
                                              deletion_t{1}},
                                              coverage_t{0, 1, 1, 0, 0, 1, 1}};
    delta_event_shared_t node6{delta_event_t{position_t{.offset = 10u},
                                              deletion_t{1}},
                                              coverage_t{0, 1, 1, 0, 0, 1, 0}};
    delta_event_shared_t node7{delta_event_t{position_t{.offset = 10u},
                                              substitution_t{"abc"s}},
                                              coverage_t{0, 1, 1, 0, 0, 1, 1}};

    EXPECT_FALSE(std::as_const(node1) < std::as_const(node1));
    EXPECT_FALSE(std::as_const(node1) < std::as_const(node2));
    EXPECT_FALSE(std::as_const(node1) < std::as_const(node3));
    EXPECT_FALSE(std::as_const(node1) < std::as_const(node4));
    EXPECT_TRUE (std::as_const(node1) < std::as_const(node5));
    EXPECT_TRUE (std::as_const(node1) < std::as_const(node6));
    EXPECT_TRUE (std::as_const(node1) < std::as_const(node7));

    EXPECT_FALSE(std::as_const(node2) < std::as_const(node1));
    EXPECT_FALSE(std::as_const(node2) < std::as_const(node2));
    EXPECT_FALSE(std::as_const(node2) < std::as_const(node3));
    EXPECT_FALSE(std::as_const(node2) < std::as_const(node4));
    EXPECT_TRUE (std::as_const(node2) < std::as_const(node5));
    EXPECT_TRUE (std::as_const(node2) < std::as_const(node6));
    EXPECT_TRUE (std::as_const(node2) < std::as_const(node7));

    EXPECT_FALSE(std::as_const(node3) < std::as_const(node1));
    EXPECT_FALSE(std::as_const(node3) < std::as_const(node2));
    EXPECT_FALSE(std::as_const(node3) < std::as_const(node3));
    EXPECT_FALSE(std::as_const(node3) < std::as_const(node4));
    EXPECT_TRUE (std::as_const(node3) < std::as_const(node5));
    EXPECT_TRUE (std::as_const(node3) < std::as_const(node6));
    EXPECT_TRUE (std::as_const(node3) < std::as_const(node7));

    EXPECT_FALSE(std::as_const(node4) < std::as_const(node1));
    EXPECT_FALSE(std::as_const(node4) < std::as_const(node2));
    EXPECT_FALSE(std::as_const(node4) < std::as_const(node3));
    EXPECT_FALSE(std::as_const(node4) < std::as_const(node4));
    EXPECT_TRUE (std::as_const(node4) < std::as_const(node5));
    EXPECT_TRUE (std::as_const(node4) < std::as_const(node6));
    EXPECT_TRUE (std::as_const(node4) < std::as_const(node7));

    EXPECT_FALSE(std::as_const(node5) < std::as_const(node1));
    EXPECT_FALSE(std::as_const(node5) < std::as_const(node2));
    EXPECT_FALSE(std::as_const(node5) < std::as_const(node3));
    EXPECT_FALSE(std::as_const(node5) < std::as_const(node4));
    EXPECT_FALSE(std::as_const(node5) < std::as_const(node5));
    EXPECT_FALSE(std::as_const(node5) < std::as_const(node6));
    EXPECT_FALSE(std::as_const(node5) < std::as_const(node7));

    EXPECT_FALSE(std::as_const(node6) < std::as_const(node1));
    EXPECT_FALSE(std::as_const(node6) < std::as_const(node2));
    EXPECT_FALSE(std::as_const(node6) < std::as_const(node3));
    EXPECT_FALSE(std::as_const(node6) < std::as_const(node4));
    EXPECT_FALSE(std::as_const(node6) < std::as_const(node5));
    EXPECT_FALSE(std::as_const(node6) < std::as_const(node6));
    EXPECT_FALSE(std::as_const(node6) < std::as_const(node7));

    EXPECT_FALSE(std::as_const(node7) < std::as_const(node2));
    EXPECT_FALSE(std::as_const(node7) < std::as_const(node3));
    EXPECT_FALSE(std::as_const(node7) < std::as_const(node1));
    EXPECT_FALSE(std::as_const(node7) < std::as_const(node4));
    EXPECT_FALSE(std::as_const(node7) < std::as_const(node5));
    EXPECT_FALSE(std::as_const(node7) < std::as_const(node6));
    EXPECT_FALSE(std::as_const(node7) < std::as_const(node7));
}

TEST_F(delta_event_shared_fixture, fill_vector)
{
    std::vector nodes
    {
        delta_event_shared_t{delta_event_t{position_t{.offset = 10u}, deletion_t{1}}, coverage_t{1, 0}},
        delta_event_shared_t{delta_event_t{position_t{.offset = 11u}, deletion_t{1}}, coverage_t{0, 1}},
    };

    nodes.emplace_back(position_t{.offset = 12u}, insertion_t{"xxx"}, coverage_t{1, 1});

    EXPECT_EQ(nodes.size(), 3u);
}

TEST_F(delta_event_shared_fixture, formatted_output)
{
    using namespace std::literals;

    delta_event_shared_t del{position_t{.offset = 10u}, deletion_t{1}, coverage_t{1, 0, 1, 0}};
    delta_event_shared_t ins{position_t{.offset = 11u}, insertion_t{"ii"s}, coverage_t{1, 1, 0, 0}};
    delta_event_shared_t sub{position_t{.offset = 12u}, substitution_t{"sss"s}, coverage_t{0, 0, 1, 1}};

    std::stringstream sstream;

    sstream << del << '\n' << ins << '\n' << sub;
    EXPECT_EQ(sstream.str(),
             ("([idx: 0, pos: 10], del: 1) ~ <1010>\n([idx: 0, pos: 11], ins: ii) ~ <1100>\n([idx: 0, pos: 12], sub: sss) ~ <0011>"s));
}

TEST_F(delta_event_shared_fixture, save_substitution)
{
    using namespace std::literals;
    std::stringstream archive_stream{};

    delta_event_shared_t substitution_event{position_t{.offset = 23u}, substitution_t{"abcd"s}, coverage_t{1, 0, 1, 0}};

    {
        cereal::JSONOutputArchive output_archive(archive_stream);
        substitution_event.save(output_archive);
    }

    EXPECT_EQ(archive_stream.str(), expected_substitution_archive);
}

TEST_F(delta_event_shared_fixture, save_insertion)
{
    using namespace std::literals;
    std::stringstream archive_stream{};

    delta_event_shared_t insertion_event{position_t{.offset = 5u}, insertion_t{"ijklm"s}, coverage_t{1, 0, 1, 0}};

    {
        cereal::JSONOutputArchive output_archive(archive_stream);
        insertion_event.save(output_archive);
    }

    EXPECT_EQ(archive_stream.str(), expected_insertion_archive);
}

TEST_F(delta_event_shared_fixture, save_deletion)
{
    std::stringstream archive_stream{};

    delta_event_shared_t deletion_event{position_t{.offset = 100u}, deletion_t{10}, coverage_t{1, 0, 1, 0}};

    {
        cereal::JSONOutputArchive output_archive(archive_stream);
        deletion_event.save(output_archive);
    }

    EXPECT_EQ(archive_stream.str(), expected_deletion_archive);
}

TEST_F(delta_event_shared_fixture, load_substitution)
{
    using namespace std::literals;

    std::stringstream archive_stream{expected_substitution_archive.data()};

    delta_event_shared_t substitution_event{};
    {
        cereal::JSONInputArchive input_archive(archive_stream);
        substitution_event.load(input_archive);
    }

    EXPECT_EQ(substitution_event,
              (delta_event_shared_t{position_t{.offset = 23u},
                                    libjst::detail::delta_kind_substitution{"abcd"s},
                                    coverage_t{1, 0, 1, 0}}));
}

TEST_F(delta_event_shared_fixture, load_insertion)
{
    using namespace std::literals;

    std::stringstream archive_stream{expected_insertion_archive.data()};

    delta_event_shared_t insertion_event{};
    {
        cereal::JSONInputArchive input_archive(archive_stream);
        insertion_event.load(input_archive);
    }

    EXPECT_EQ(insertion_event,
              (delta_event_shared_t{position_t{.offset = 5u},
                                    libjst::detail::delta_kind_insertion{"ijklm"s},
                                    coverage_t{1, 0, 1, 0}}));
}

TEST_F(delta_event_shared_fixture, load_deletion)
{
    std::stringstream archive_stream{expected_deletion_archive.data()};

    delta_event_shared_t deletion_event{};
    {
        cereal::JSONInputArchive input_archive(archive_stream);
        deletion_event.load(input_archive);
    }

    EXPECT_EQ(deletion_event,
              (delta_event_shared_t{position_t{.offset = 100u},
                                    libjst::detail::delta_kind_deletion{10},
                                    coverage_t{1, 0, 1, 0}}));
}
