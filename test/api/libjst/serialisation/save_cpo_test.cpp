// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

// add include for std::as_const
#include <utility>

#include <libjst/serialisation/concept.hpp>

#if __has_include(<cereal/cereal.hpp>)
#include <cereal/archives/binary.hpp>
#endif

namespace libjst::test
{

    struct test_object_free_save
    {
        int a{};
        int b{};
        int c{};
    };

    template <typename archive_t>
    void save(test_object_free_save const &object, archive_t &archive)
    {
        archive(object.a, object.b, object.c);
    }
} // namespace libjst::test

TEST(save_cpo, using_free_function)
{
#if __has_include(<cereal/cereal.hpp>)
    {

        std::stringstream stream{};
        {
            cereal::BinaryOutputArchive oarchive{stream};
            libjst::test::test_object_free_save object{.a = 1, .b = 2, .c = 3};
            libjst::save(std::as_const(object), oarchive);
        }

        libjst::test::test_object_free_save expected_object{};
        {
            cereal::BinaryInputArchive iarchive{stream};
            iarchive(expected_object.a, expected_object.b, expected_object.c);
        }

        EXPECT_EQ(expected_object.a, 1);
        EXPECT_EQ(expected_object.b, 2);
        EXPECT_EQ(expected_object.c, 3);
    }
#else
    {
        // No cereal available -> save should throw
        std::stringstream archive_mock{};
        libjst::test::test_object_free_save object{.a = 1, .b = 2, .c = 3};
        EXPECT_THROW(libjst::save(std::as_const(object), archive_mock), std::runtime_error);
    }
#endif
}

namespace libjst::test
{
    struct test_object_member_save
    {
        int a{};
        int b{};
        int c{};

        template <typename archive_t>
        void save(archive_t &archive) const
        {
            archive(a, b, c);
        }
    };
} // namespace libjst::test

TEST(save_cpo, using_member_function)
{
#if __has_include(<cereal/cereal.hpp>)
    {
        std::stringstream stream{};
        {
            cereal::BinaryOutputArchive oarchive{stream};
            libjst::test::test_object_member_save object{.a = 1, .b = 2, .c = 3};
            std::as_const(object).save(oarchive);
        }

        libjst::test::test_object_member_save expected_object{};
        {
            cereal::BinaryInputArchive iarchive{stream};
            iarchive(expected_object.a, expected_object.b, expected_object.c);
        }

        EXPECT_EQ(expected_object.a, 1);
        EXPECT_EQ(expected_object.b, 2);
        EXPECT_EQ(expected_object.c, 3);
    }
#else
    {
        // No cereal available -> save should throw
        std::stringstream archive_mock{};
        libjst::test::test_object_member_save object{.a = 1, .b = 2, .c = 3};
        EXPECT_THROW(std::as_const(object).save(archive_mock), std::runtime_error);
    }
#endif
}

namespace libjst::test
{
    struct test_object_save_tag_invoke
    {
        int a{};
        int b{};
        int c{};

    private:
        template <typename archive_t>
        friend void tag_invoke(std::tag_t<libjst::save>, test_object_save_tag_invoke const &object, archive_t &archive)
        {
            archive(object.a, object.b, object.c);
        }
    };

} // namespace libjst::test

TEST(save_cpo, using_tag_invoke)
{
#if __has_include(<cereal/cereal.hpp>)
    {
        std::stringstream stream{};
        {
            cereal::BinaryOutputArchive oarchive{stream};
            libjst::test::test_object_save_tag_invoke object{.a = 1, .b = 2, .c = 3};
            libjst::save(std::as_const(object), oarchive);
        }

        libjst::test::test_object_save_tag_invoke expected_object{};
        {
            cereal::BinaryInputArchive iarchive{stream};
            iarchive(expected_object.a, expected_object.b, expected_object.c);
        }

        EXPECT_EQ(expected_object.a, 1);
        EXPECT_EQ(expected_object.b, 2);
        EXPECT_EQ(expected_object.c, 3);
    }
#else
    {
        // No cereal available -> save should throw
        std::stringstream archive_mock{};
        libjst::test::test_object_save_tag_invoke object{.a = 1, .b = 2, .c = 3};
        EXPECT_THROW(libjst::save(std::as_const(object), archive_mock), std::runtime_error);
    }
#endif
}
