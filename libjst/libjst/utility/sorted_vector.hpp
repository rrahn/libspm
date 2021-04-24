// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::sorted_vector.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <compare>
#include <concepts>
#include <memory_resource>
#include <vector>

namespace libjst
{
template <std::semiregular key_t, typename compare_t = std::less<key_t>>
class sorted_vector
{
    using container_t = std::pmr::vector<key_t>;

    container_t _elements; //!< The container holding the elements.
public:

    template <bool is_const>
    class bi_iterator;

    using value_type = std::ranges::range_value_t<container_t>;
    using iterator = bi_iterator<false>;
    using const_iterator = bi_iterator<true>;
    using size_type = size_t;

    /*!\name Iterators
     * \{
     */
    constexpr iterator begin() noexcept
    {
        return iterator{_elements.begin()};
    }

    constexpr const_iterator begin() const noexcept
    {
        return const_iterator{_elements.begin()};
    }

    constexpr iterator end() noexcept
    {
        return iterator{_elements.end()};
    }

    constexpr const_iterator end() const noexcept
    {
        return const_iterator{_elements.end()};
    }
    //!\}

    /*!\name Cpacity
     * \{
     */
    constexpr bool empty() const noexcept
    {
        return _elements.empty();
    }

    constexpr size_type size() const noexcept
    {
        return _elements.size();
    }

    constexpr size_type max_size() const noexcept
    {
        return _elements.max_size();
    }
    //!\}

    /*!\name Modifiers
     * \{
     */
    constexpr auto clear() noexcept
    {
        _elements.clear();
    }

    iterator insert(value_type const & value)
    {
        return insert(end(), value);
    }

    iterator insert(value_type && value)
    {
        return insert(end(), std::move(value));
    }

    iterator insert(const_iterator hint, value_type const & value)
    {
        return insert_impl(hint, value);
    }

    iterator insert(const_iterator hint, value_type && value)
    {
        return insert_impl(hint, std::move(value));
    }

    template <typename ...args_t>
    iterator emplace(args_t && ...args)
    {
        return insert(value_type(std::forward<args_t>(args)...));
    }

    template <typename ...args_t>
    iterator emplace_hint(const_iterator hint, args_t && ...args)
    {
        return insert(hint, value_type(std::forward<args_t>(args)...));
    }

    iterator erase(iterator pos)
    {
        return iterator{_elements.erase(pos.base())};
    }

    iterator erase(const_iterator pos)
    {
        return iterator{_elements.erase(pos.base())};
    }

    iterator erase(const_iterator first, const_iterator last)
    {
        return iterator{_elements.erase(first.base(), last.base())};
    }

    size_type erase(key_t const & key)
    {
        auto [first, last] = equal_range_impl(key);
        size_t erased_elements = std::ranges::distance(first, last);
        _elements.erase(first.base(), last.base());
        return erased_elements;
    }
    // constexpr void swap();
    // constexpr void extract();
    // constexpr void merge();
    //!\}

    /*!\name Lookup
     * \{
     */
    size_type count(key_t const & key) const
    {
        auto [first, last] = equal_range(key);
        return std::ranges::distance(first, last);
    }

    template<typename other_key_t>
    size_type count(other_key_t const & key) const
    {
        auto [first, last] = equal_range(key);
        return std::ranges::distance(first, last);
    }

    iterator find(key_t const & key)
    {
        return find_impl(key);
    }

    const_iterator find(key_t const & key) const
    {
        return find_impl(key);
    }

    template<typename comparable_key_t>
    iterator find(comparable_key_t const & key)
    {
        return find_impl(key);
    }

    template<typename comparable_key_t>
    const_iterator find(comparable_key_t const & key) const
    {
        return find_impl(key);
    }

    bool contains(key_t const & key) const
    {
        return find(key) != end();
    }

    template<typename comparable_key_t>
    bool contains(comparable_key_t const & key) const
    {
        return find(key) != end();
    }

    std::pair<iterator, iterator> equal_range(key_t const & key)
    {
        return equal_range_impl(key);
    }
    std::pair<const_iterator, const_iterator> equal_range(key_t const & key) const
    {
        return equal_range_impl(key);
    }

    template<typename other_key_t>
    std::pair<iterator, iterator> equal_range(other_key_t const & key)
    {
        return equal_range_impl(key);
    }

    template<typename other_key_t>
    std::pair<const_iterator, const_iterator> equal_range(other_key_t const & key) const
    {
        return equal_range_impl(key);
    }

    // iterator lower_bound(key_t const & key)
    // {
    //     return lower_bound_impl(key);
    // }

    // const_iterator lower_bound(key_t const & key) const
    // {
    //     return lower_bound_impl(key);
    // }

    template<typename other_key_t>
        requires std::strict_weak_order<compare_t, key_t, other_key_t>
    iterator lower_bound(other_key_t const & key)
    {
        return lower_bound_impl(key);
    }

    template<typename other_key_t>
        requires std::strict_weak_order<compare_t, key_t, other_key_t>
    const_iterator lower_bound(other_key_t const & key) const
    {
        return lower_bound_impl(key);
    }

    iterator lower_bound(key_t const & key)
    {
        return lower_bound_impl(key);
    }
    const_iterator lower_bound(key_t const & key) const
    {
        return lower_bound_impl(key);
    }

    template<typename other_key_t>
    iterator upper_bound(other_key_t const & key)
    {
        return upper_bound_impl(key);
    }

    template<typename other_key_t>
    const_iterator upper_bound(other_key_t const & key) const
    {
        return upper_bound_impl(key);
    }
    //!\}

    /*!\name Comparison
     * \{
     */
    bool operator==(sorted_vector const & other) const = default;
    std::strong_ordering operator<=>(sorted_vector const & other) const = default;
    //!\}

private:

    template <typename comparable_key_t>
    iterator find_impl(comparable_key_t && key)
    {
        auto base_it = std::ranges::lower_bound(_elements, key, compare_t{});
        if (compare_t{}(key, *base_it)) // not identical
            return iterator{_elements.end()};

        return iterator{base_it};
    }

    template <typename comparable_key_t>
    const_iterator find_impl(comparable_key_t && key) const
    {
        auto base_it = std::ranges::lower_bound(_elements, key, compare_t{});
        if (compare_t{}(key, *base_it)) // not identical
            return const_iterator{_elements.end()};

        return const_iterator{base_it};
    }

    template <typename value_t>
    iterator insert_impl(const_iterator hint, value_t && value)
    {
        // We want to insert before the hint but only if it is ok to insert here.
        if (_elements.empty())
        {
            _elements.push_back(std::forward<value_t>(value));
            return iterator{_elements.begin()};
        }

        compare_t compare{};

        auto hint_base =  hint.base();
        bool const at_end = hint_base == _elements.end();
        // (at_end or value < value at hint) && value >= value directly before hint
        hint_base = ((at_end || compare(value, *hint_base)) && (hint_base == _elements.begin() || !compare(value, *(std::ranges::prev(hint_base)))))
                    ? hint_base
                    : std::ranges::upper_bound(_elements, value, compare);

        return iterator{_elements.insert(hint_base, std::forward<value_t>(value))};
    }

    template <typename comparable_key_t>
    std::pair<iterator, iterator> equal_range_impl(comparable_key_t const & key)
    {
        auto rng = std::ranges::equal_range(_elements, key, compare_t{});
        return {iterator{rng.begin()}, iterator{rng.end()}};
    }

    template <typename comparable_key_t>
    std::pair<const_iterator, const_iterator> equal_range_impl(comparable_key_t const & key) const
    {
        auto rng = std::ranges::equal_range(_elements, key, compare_t{});
        return {const_iterator{rng.begin()}, const_iterator{rng.end()}};
    }

    template <typename comparable_key_t>
    iterator lower_bound_impl(comparable_key_t const & key)
    {
        return iterator{std::ranges::lower_bound(_elements, key, compare_t{})};
    }

    template <typename comparable_key_t>
    const_iterator lower_bound_impl(comparable_key_t const & key) const
    {
        return const_iterator{std::ranges::lower_bound(_elements, key, compare_t{})};
    }

    template <typename comparable_key_t>
    iterator upper_bound_impl(comparable_key_t const & key)
    {
        return iterator{std::ranges::upper_bound(_elements, key, compare_t{})};
    }

    template <typename comparable_key_t>
    const_iterator upper_bound_impl(comparable_key_t const & key) const
    {
        return const_iterator{std::ranges::upper_bound(_elements, key, compare_t{})};
    }
};

template <std::semiregular key_t, typename compare_t>
template <bool is_const>
class sorted_vector<key_t, compare_t>::bi_iterator
{
private:
    using maybe_const_container_type = std::conditional_t<is_const, container_t const, container_t>;
    using iterator_type = std::ranges::iterator_t<maybe_const_container_type>;

    template <bool>
    friend class bi_iterator;

    iterator_type _iter;
public:

    using value_type = std::ranges::range_value_t<maybe_const_container_type>;
    using reference = std::ranges::range_reference_t<maybe_const_container_type>;
    using difference_type = std::ranges::range_difference_t<maybe_const_container_type>;
    using pointer = std::conditional_t<is_const, typename container_t::const_pointer, typename container_t::pointer>;
    using iterator_category = std::bidirectional_iterator_tag;

    constexpr bi_iterator() = default;
    constexpr explicit bi_iterator(iterator_type iter) noexcept : _iter{std::move(iter)}
    {}

    constexpr bi_iterator(bi_iterator<!is_const> iter) noexcept
        requires is_const
        : _iter{std::move(iter._iter)}
    {}

    iterator_type base() const &
    {
        return _iter;
    }

    iterator_type base() &&
    {
        return std::move(_iter);
    }

    constexpr reference operator*() const noexcept
    {
        return *_iter;
    }

    constexpr pointer operator->() const noexcept
    {
        return _iter.operator->();
    }

    constexpr bi_iterator & operator++() noexcept
    {
        ++_iter;
        return *this;
    }

    constexpr bi_iterator operator++(int) noexcept
    {
        iterator tmp{*this};
        ++(*this);
        return tmp;
    }

    constexpr bi_iterator & operator--() noexcept
    {
        --_iter;
        return *this;
    }

    constexpr bi_iterator operator--(int) noexcept
    {
        iterator tmp{*this};
        --(*this);
        return tmp;
    }

    template <bool other_const>
    constexpr difference_type operator-(bi_iterator<other_const> const & rhs) const noexcept
    {
        return _iter - rhs._iter;
    }

    template <bool other_const>
    constexpr bool operator==(bi_iterator<other_const> const & rhs) const noexcept
    {
        return _iter == rhs._iter;
    }

    template <bool other_const>
    constexpr std::strong_ordering operator<=>(bi_iterator<other_const> const & rhs) const noexcept
    {
        return _iter <=> rhs._iter;
    }
};

}  // namespace libjst
