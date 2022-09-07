// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <concepts>
#include <ranges>
#include <type_traits>

#include <seqan/basic.h>

namespace seqan
{
    namespace detail
    {
        template <typename range_t>
        concept adaptable_view = Not<Is<ContainerConcept<std::remove_reference_t<range_t> > > >::VALUE &&
                                 std::ranges::view<range_t> &&
                                 std::ranges::common_range<range_t>;

        template <typename range_t>
        struct range_traits
        {
            using value_type = std::ranges::range_value_t<range_t>;
            using reference = std::ranges::range_reference_t<range_t>;
            using difference_type = std::ranges::range_difference_t<range_t>;
        };
    } // namespace detail

    // ----------------------------------------------------------------------------
    // Overloading StdContainerIterator
    // ----------------------------------------------------------------------------

    template <typename range_t>
        requires detail::adaptable_view<range_t>
    struct StdContainerIterator<range_t>
    {
        using Type = std::ranges::iterator_t<range_t>;
    };

    template <typename range_t>
        requires detail::adaptable_view<range_t>
    struct StdContainerIterator<range_t const>
    {
        using Type = std::ranges::iterator_t<range_t const>;
    };

    // ----------------------------------------------------------------------------
    // Overloading Iterator
    // ----------------------------------------------------------------------------

    // standard
    template <typename range_t>
        requires detail::adaptable_view<range_t>
    struct Iterator<range_t, Standard>
    {
    private:
        using range_object_t = std::remove_reference_t<range_t>;
    public:
        using Type = Iter<range_object_t, StdIteratorAdaptor>;
    };

    // standard & const
    template <typename range_t>
        requires detail::adaptable_view<range_t const>
    struct Iterator<range_t const, Standard>
    {
    private:
        using range_object_t = std::remove_reference_t<range_t const>;
    public:
        using Type = Iter<range_object_t, StdIteratorAdaptor>;
    };

    // rooted
    template <typename range_t>
        requires detail::adaptable_view<range_t>
    struct Iterator<range_t, Rooted>
    {
    private:
        using range_object_t = std::remove_reference_t<range_t>;
    public:
        using Type = Iter<range_object_t, AdaptorIterator<Iter<range_object_t, StdIteratorAdaptor> > >;
    };

    // rooted & const
    template <typename range_t>
        requires detail::adaptable_view<range_t const>
    struct Iterator<range_t const, Rooted>
    {
    private:
        using range_object_t = std::remove_reference_t<range_t const>;
    public:
        using Type = Iter<range_object_t, AdaptorIterator<Iter<range_object_t, StdIteratorAdaptor> > >;
    };

    // ----------------------------------------------------------------------------
    // Overloading range type traits
    // ----------------------------------------------------------------------------

    // // Value
    // template <typename range_t>
    //     requires detail::adaptable_view<range_t>
    // struct Value<range_t>
    // {
    //     using Type = typename detail::range_traits<range_t>::value_type;
    // };

    // // Value & const
    // template <typename range_t>
    //     requires detail::adaptable_view<range_t const>
    // struct Value<range_t const>
    // {
    //     using Type = typename detail::range_traits<range_t const>::value_type;
    // };

    // // Reference
    template <typename range_t>
        requires detail::adaptable_view<range_t>
    struct Reference<range_t>
    {
        using Type = typename detail::range_traits<range_t>::reference;
    };

    // Reference & const
    template <typename range_t>
        requires detail::adaptable_view<range_t const>
    struct Reference<range_t const>
    {
        using Type = typename detail::range_traits<range_t const>::reference;
    };

    // // Position
    // template <typename range_t>
    //     requires detail::adaptable_view<range_t>
    // struct Position<range_t>
    // {
    //     using Type = typename detail::range_traits<range_t>::position_type;
    // };

    // // Position & const
    // template <typename range_t>
    //     requires detail::adaptable_view<range_t const>
    // struct Position<range_t const>
    // {
    //     using Type = typename detail::range_traits<range_t const>::position_type;
    // };

    // Value
    template <typename range_t>
        requires detail::adaptable_view<range_t>
    struct Value<Iter<range_t, AdaptorIterator<Iter<range_t, StdIteratorAdaptor> > > > :
        public Value<Iter<range_t, StdIteratorAdaptor>>
    {
    };

    // Value & const
    template <typename range_t>
        requires detail::adaptable_view<range_t const>
    struct Value<Iter<range_t const, AdaptorIterator<Iter<range_t const, StdIteratorAdaptor> > > > :
        public Value<Iter<range_t const, StdIteratorAdaptor>>
    {
    };

    // Reference
    // inline typename Reference<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > >::Type
    template <typename range_t, typename spec_t>
        requires detail::adaptable_view<range_t>
    struct Reference<Iter<range_t, AdaptorIterator<Iter<range_t, StdIteratorAdaptor>, spec_t> > > :
        public Reference<Iter<range_t, StdIteratorAdaptor>>
    {
    };

    // Reference & const
    template <typename range_t, typename spec_t>
        requires detail::adaptable_view<range_t const>
    struct Reference<Iter<range_t const, AdaptorIterator<Iter<range_t const, StdIteratorAdaptor>, spec_t> > > :
        public Reference<Iter<range_t const, StdIteratorAdaptor>>
    {
    };

    // Position
    template <typename range_t>
        requires detail::adaptable_view<range_t>
    struct Position<Iter<range_t, AdaptorIterator<Iter<range_t, StdIteratorAdaptor> > > > :
        public Position<Iter<range_t, StdIteratorAdaptor>>
    {
    };

    // Position & const
    template <typename range_t>
        requires detail::adaptable_view<range_t const>
    struct Position<Iter<range_t const, AdaptorIterator<Iter<range_t const, StdIteratorAdaptor> > > > :
        public Position<Iter<range_t const, StdIteratorAdaptor>>
    {
    };

    // ----------------------------------------------------------------------------
    // Overloading begin
    // ----------------------------------------------------------------------------

    // Standard
    template <typename range_t>
        requires detail::adaptable_view<range_t>
    inline auto begin(range_t &range, Standard)
        -> typename Iterator<range_t, Standard>::Type
    {
        using iter_t = typename Iterator<range_t, Standard>::Type;
        return iter_t{std::ranges::begin(range)};
    }

    // Standard & const
    template <typename range_t>
        requires detail::adaptable_view<range_t>
    inline auto begin(range_t const &range, Standard)
        -> typename Iterator<range_t const, Standard>::Type
    {
        using iter_t = typename Iterator<range_t const, Standard>::Type;
        return iter_t{std::ranges::begin(range)};
    }

    // Rooted
    template <typename range_t>
        requires detail::adaptable_view<range_t>
    inline auto begin(range_t &range, Rooted)
        -> typename Iterator<range_t, Rooted>::Type
    {
        std::cout << "Expected begin\n";
        using iter_t = typename Iterator<range_t, Rooted>::Type;
        static_assert(std::same_as<iter_t, Iter<range_t, AdaptorIterator<Iter<range_t, StdIteratorAdaptor>>>>, "wrong begin iterator");
        return iter_t{range, begin(range, Standard{})};
    }

    // Rooted & const
    template <typename range_t>
        requires detail::adaptable_view<range_t>
    inline auto begin(range_t const &range, Rooted)
        -> typename Iterator<range_t const, Rooted>::Type
    {
        using iter_t = typename Iterator<range_t const, Rooted>::Type;
        return iter_t{range, begin(range, Standard{})};
    }

    // ----------------------------------------------------------------------------
    // Overloading end
    // ----------------------------------------------------------------------------

    // Standard
    template <typename range_t>
        requires detail::adaptable_view<range_t>
    inline auto end(range_t &range, Standard)
        -> typename Iterator<range_t, Standard>::Type
    {
        using iter_t = typename Iterator<range_t, Standard>::Type;
        return iter_t{std::ranges::end(range)};
    }

    // Standard & const
    template <typename range_t>
        requires detail::adaptable_view<range_t>
    inline auto end(range_t const &range, Standard)
        -> typename Iterator<range_t const, Standard>::Type
    {
        using iter_t = typename Iterator<range_t const, Standard>::Type;
        return iter_t{std::ranges::end(range)};
    }

    // Rooted
    template <typename range_t>
        requires detail::adaptable_view<range_t>
    inline auto end(range_t &range, Rooted)
        -> typename Iterator<range_t, Rooted>::Type
    {
        std::cout << "Expected end\n";
        using iter_t = typename Iterator<range_t, Rooted>::Type;
        static_assert(std::same_as<iter_t, Iter<range_t, AdaptorIterator<Iter<range_t, StdIteratorAdaptor>>>>, "wrong end iterator");
        return iter_t{range, end(range, Standard{})};
    }

    // Rooted & const
    template <typename range_t>
        requires detail::adaptable_view<range_t>
    inline auto end(range_t const &range, Rooted)
        -> typename Iterator<range_t const, Rooted>::Type
    {
        using iter_t = typename Iterator<range_t const, Rooted>::Type;
        return iter_t{range, end(range, Standard{})};
    }

} // namespace seqan
