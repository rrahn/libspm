// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides generic random access iterator for variant stores.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <memory>

namespace libjst
{

template <typename variant_store_t>
class variant_store_iterator
{
private:
    variant_store_t const *_variant_store{};
    uint64_t _position{};

public:
    using value_type = typename variant_store_t::value_type;
    using reference = typename variant_store_t::reference;
    using difference_type =	std::ptrdiff_t;
    using pointer =	void;
    using iterator_category	= std::random_access_iterator_tag;

    constexpr variant_store_iterator() = default;
    constexpr explicit variant_store_iterator(variant_store_t const & variant_store, size_t const position) noexcept :
        _variant_store{std::addressof(variant_store)},
        _position{position}
    {}

    constexpr reference operator*() const noexcept
    {
        return (*_variant_store)[_position];
    }

    constexpr reference operator[](difference_type const offset) const noexcept
    {
        return *(*this + offset);
    }

    constexpr variant_store_iterator & operator++() noexcept
    {
        ++_position;
        return *this;
    }

    constexpr variant_store_iterator operator++(int)
        noexcept(std::is_nothrow_copy_constructible_v<variant_store_iterator>)
    {
        variant_store_iterator tmp{*this};
        ++(*this);
        return tmp;
    }

    constexpr variant_store_iterator & operator+=(difference_type const offset) noexcept
    {
        _position += offset;
        return *this;
    }

    constexpr variant_store_iterator operator+(difference_type const offset) const
        noexcept(std::is_nothrow_copy_constructible_v<variant_store_iterator>)
    {
        variant_store_iterator tmp{*this};
        return tmp += offset;
    }

    friend constexpr variant_store_iterator operator+(difference_type const offset, variant_store_iterator const & rhs)
        noexcept(noexcept(rhs + offset))
    {
        return rhs + offset;
    }

    constexpr variant_store_iterator & operator--() noexcept
    {
        --_position;
        return *this;
    }

    constexpr variant_store_iterator operator--(int)
        noexcept(std::is_nothrow_copy_constructible_v<variant_store_iterator>)
    {
        variant_store_iterator tmp{*this};
        --(*this);
        return tmp;
    }

    constexpr variant_store_iterator & operator-=(difference_type const offset) noexcept
    {
        _position -= offset;
        return *this;
    }

    constexpr variant_store_iterator operator-(difference_type const offset) const
        noexcept(std::is_nothrow_copy_constructible_v<variant_store_iterator>)
    {
        variant_store_iterator tmp{*this};
        return tmp -= offset;
    }

    constexpr difference_type operator-(variant_store_iterator const & rhs) const noexcept
    {
        return _position - rhs._position;
    }

    constexpr bool operator==(variant_store_iterator const &) const noexcept = default;

    constexpr std::strong_ordering operator<=>(variant_store_iterator const & rhs) const noexcept
    {
        return _position <=> rhs._position;
    }
};
}  // namespace libjst
