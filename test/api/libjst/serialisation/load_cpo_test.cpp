//// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <libjst/serialisation/concept.hpp>

// Check if cereal is available
#if __has_include(<cereal/cereal.hpp>)
#include <cereal/archives/binary.hpp>
#endif

namespace libjst::test
{

    struct test_object_free_load
    {
        int a{};
        int b{};
        int c{};
    };

    template <typename archive_t>
    void load(test_object_free_load &object, archive_t &archive)
    {
        archive(object.a, object.b, object.c);
    }
} // namespace libjst::test

TEST(load_cpo, using_free_function)
{
#if __has_include(<cereal/cereal.hpp>)
    {
        std::stringstream stream{};
        {
            cereal::BinaryOutputArchive oarchive{stream};
            oarchive(1, 2, 3);
        }

        libjst::test::test_object_free_load object{};
        // Test before loading
        EXPECT_EQ(object.a, 0);
        EXPECT_EQ(object.b, 0);
        EXPECT_EQ(object.c, 0);
        {
            cereal::BinaryInputArchive iarchive{stream};
            libjst::load(object, iarchive);
        }
        // Test after loading
        EXPECT_EQ(object.a, 1);
        EXPECT_EQ(object.b, 2);
        EXPECT_EQ(object.c, 3);
    }
#else
    {
        libjst::test::test_object_free_load test_object{};
        std::stringstream archive_mock{"1, 2, 3"};
        EXPECT_THROW(libjst::load(test_object, archive_mock), std::runtime_error);
    }
#endif
}

namespace libjst::test
{

    struct test_object_member_load
    {
        int a{};
        int b{};
        int c{};

        template <typename archive_t>
        void load(archive_t &archive)
        {
            archive(a, b, c);
        }
    };
}

TEST(load_cpo, using_member_function)
{
#if __has_include(<cereal/cereal.hpp>)
    {
        std::stringstream stream{};
        {
            cereal::BinaryOutputArchive oarchive{stream};
            oarchive(1, 2, 3);
        }

        libjst::test::test_object_member_load object{};
        // Test before loading
        EXPECT_EQ(object.a, 0);
        EXPECT_EQ(object.b, 0);
        EXPECT_EQ(object.c, 0);
        {
            cereal::BinaryInputArchive iarchive{stream};
            libjst::load(object, iarchive);
        }
        // Test after loading
        EXPECT_EQ(object.a, 1);
        EXPECT_EQ(object.b, 2);
        EXPECT_EQ(object.c, 3);
    }
#else
    {
        libjst::test::test_object_member_load test_object{};
        std::stringstream archive_mock{"1, 2, 3"};
        EXPECT_THROW(libjst::load(test_object, archive_mock), std::runtime_error);
    }
#endif
}

namespace libjst::test
{

    struct test_object_tag_invoke
    {
        int a{};
        int b{};
        int c{};

    private:

        template <typename archive_t>
        friend void tag_invoke(std::tag_t<libjst::load>, test_object_tag_invoke &object, archive_t &archive)
        {
            archive(object.a, object.b, object.c);
        }
    };
} // namespace libjst::test

TEST(load_cpo, using_tag_invoke_friend)
{
#if __has_include(<cereal/cereal.hpp>)
    {
        std::stringstream stream{};
        {
            cereal::BinaryOutputArchive oarchive{stream};
            oarchive(1, 2, 3);
        }

        libjst::test::test_object_tag_invoke object{};
        // Test before loading
        EXPECT_EQ(object.a, 0);
        EXPECT_EQ(object.b, 0);
        EXPECT_EQ(object.c, 0);
        {
            cereal::BinaryInputArchive iarchive{stream};
            libjst::load(object, iarchive);
        }
        // Test after loading
        EXPECT_EQ(object.a, 1);
        EXPECT_EQ(object.b, 2);
        EXPECT_EQ(object.c, 3);
    }
#else
    {
        libjst::test::test_object_tag_invoke test_object{};
        std::stringstream archive_mock{"1, 2, 3"};
        EXPECT_THROW(libjst::load(test_object, archive_mock), std::runtime_error);
    }
#endif
}
