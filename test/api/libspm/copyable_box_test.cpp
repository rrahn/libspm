// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <libspm/copyable_box.hpp>

namespace test
{
    template <typename t, bool has_default_init, bool is_copyable, bool is_movable>
    class non_default_type
    {
    private:
        t _value;
    public:

        constexpr non_default_type()
            requires has_default_init
            : _value{}
        {}
        constexpr non_default_type(non_default_type const &) = default;
        constexpr non_default_type(non_default_type && ) = default;
        constexpr non_default_type & operator=(non_default_type const & other)
            requires is_copyable
        {
            _value = other._value;
            return *this;
        }

        constexpr non_default_type & operator=(non_default_type && other)
            requires is_movable
        {
            _value = std::move(other._value);
            return *this;
        }

        constexpr non_default_type(t value) : _value{value}
        {}
        constexpr non_default_type & operator=(t value)
        {
            _value = value;
            return *this;
        }

        constexpr operator t() const
        {
            return _value;
        }
    };


} // namespace test

template <typename t>
struct copyable_box_test : ::testing::Test
{
    using type = t;
    using box_t = spm::copyable_box<t>;
};

using test_types = ::testing::Types<int,
                                    test::non_default_type<int, false, true, true>,
                                    test::non_default_type<int, true, false, true>,
                                    test::non_default_type<int, true, true, false>,
                                    test::non_default_type<int, true, false, false>
                                >;
TYPED_TEST_SUITE(copyable_box_test, test_types);

TYPED_TEST(copyable_box_test, default_construction)
{
    EXPECT_EQ(std::default_initializable<typename TestFixture::box_t>,
              std::default_initializable<typename TestFixture::type>);
}

TYPED_TEST(copyable_box_test, copy_construction)
{
    // value construction
    typename TestFixture::box_t box{10};
    typename TestFixture::box_t box2{box};

    EXPECT_EQ(*box, *box2);
}

TYPED_TEST(copyable_box_test, move_construction)
{
    typename TestFixture::box_t box{10};
    typename TestFixture::box_t box2{std::move(box)};

    EXPECT_EQ(*box2, 10);
}

TYPED_TEST(copyable_box_test, copy_assignment)
{
    typename TestFixture::box_t box{10};
    typename TestFixture::box_t box2{1};

    EXPECT_EQ(*box2, 1);
    box2 = box;
    EXPECT_EQ(*box2, 10);
}

TYPED_TEST(copyable_box_test, move_assignment)
{
    typename TestFixture::box_t box{10};
    typename TestFixture::box_t box2{1};

    EXPECT_EQ(*box2, 1);
    box2 = std::move(box);
    EXPECT_EQ(*box2, 10);
}

TYPED_TEST(copyable_box_test, value_construction)
{
    typename TestFixture::box_t box1{10};
    typename TestFixture::box_t box2{std::in_place, 10};

    EXPECT_EQ(*box1, 10);
    EXPECT_EQ(*box2, 10);
}

TYPED_TEST(copyable_box_test, value_assignment)
{
    typename TestFixture::box_t box1{10};

    EXPECT_EQ(*box1, 10);
    box1 = 20;
    EXPECT_EQ(*box1, 20);
}

TYPED_TEST(copyable_box_test, bool_conversion)
{
    typename TestFixture::box_t box1{10};
    EXPECT_TRUE(box1);

    if constexpr (std::default_initializable<typename TestFixture::box_t>) {
        typename TestFixture::box_t box2{};
        EXPECT_TRUE(box2);
    }
}

TYPED_TEST(copyable_box_test, has_value)
{
    typename TestFixture::box_t box1{10};
    EXPECT_TRUE(box1.has_value());

    if constexpr (std::default_initializable<typename TestFixture::box_t>) {
        typename TestFixture::box_t box2{};
        EXPECT_TRUE(box2.has_value());
    }
}

TYPED_TEST(copyable_box_test, dereference)
{
    typename TestFixture::box_t box1{10};

    EXPECT_EQ(*box1, 10);
}

TYPED_TEST(copyable_box_test, reset)
{
    typename TestFixture::box_t box1{10};

    EXPECT_EQ(*box1, 10);
    box1.reset();
    EXPECT_FALSE(box1.has_value());
}

TYPED_TEST(copyable_box_test, emplace)
{
    typename TestFixture::box_t box1{1};
    EXPECT_EQ(*box1, 1);
    box1.emplace(10);
    EXPECT_EQ(*box1, 10);
}


