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
#include <ranges>
#include <vector>

#include <cereal/types/vector.hpp>

#include <seqan3/core/concept/cereal.hpp>

#include <libjst/utility/stable_random_access_iterator.hpp>

namespace libjst
{

// Forward declaration for friend access
template <std::ranges::view segment_t>
class journal_decorator;

template <typename journal_decorator_t>
class journal_decorator_revertable;

// TODO: Call this the journal
    // can be customised with the dictionary type.
    // if not specified uses its own dictionary type.
template <std::semiregular key_t, typename compare_t = std::less<key_t>>
class sorted_vector
{
    using container_t = std::vector<key_t>;

    container_t _elements; //!< The container holding the elements.

    template <bool is_const>
    class bi_iterator;

    template <std::ranges::view>
    friend class journal_decorator;

    template <typename>
    friend class journal_decorator_revertable;

public:

    using value_type = std::ranges::range_value_t<container_t>;
    using iterator = stable_random_access_iterator<container_t>;
    using const_iterator = stable_random_access_iterator<container_t const>;
    using size_type = size_t;
    using key_compare = compare_t;

    sorted_vector() noexcept = default;
    sorted_vector(sorted_vector const &) noexcept = default;
    sorted_vector(sorted_vector &&) noexcept = default;
    sorted_vector & operator=(sorted_vector const &) noexcept = default;
    sorted_vector & operator=(sorted_vector &&) noexcept = default;
    ~sorted_vector() = default;

    /*!\name Iterators
     * \{
     */
    constexpr iterator begin() noexcept
    {
        return iterator{std::addressof(_elements), 0};
    }

    constexpr const_iterator begin() const noexcept
    {
        return const_iterator{std::addressof(_elements), 0};
    }

    constexpr iterator end() noexcept
    {
        return iterator{std::addressof(_elements), std::ranges::ssize(_elements)};
    }

    constexpr const_iterator end() const noexcept
    {
        return const_iterator{std::addressof(_elements), std::ranges::ssize(_elements)};
    }
    //!\}

    constexpr container_t & data() noexcept {
        return _elements;
    }

    constexpr container_t const & data() const noexcept {
        return _elements;
    }

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

    constexpr void reserve(size_type const new_capacity)
    {
        _elements.reserve(new_capacity);
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

    // ----------------------------------------------------------------------------
    // Serialisation
    // ----------------------------------------------------------------------------

    template <seqan3::cereal_input_archive archive_t>
    void load(archive_t & iarchive)
    {
        iarchive(_elements);
    }

    template <seqan3::cereal_output_archive archive_t>
    void save(archive_t & oarchive) const
    {
        oarchive(_elements);
    }

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
            return begin();
        }

        compare_t compare{};

        auto hint_base =  hint.base();
        bool const at_end = hint_base == _elements.end();
        // (at_end or value < value at hint) && value >= value directly before hint
        hint_base = ((at_end || compare(value, *hint_base)) && (hint_base == _elements.begin() || !compare(value, *(std::ranges::prev(hint_base)))))
                    ? hint_base
                    : std::ranges::upper_bound(_elements, value, compare);

        auto insert_it = _elements.insert(hint_base, std::forward<value_t>(value));
        auto insert_pos = std::ranges::distance(_elements.begin(), insert_it);
        return iterator{std::addressof(_elements), insert_pos};
    }

    template <typename comparable_key_t>
    std::pair<iterator, iterator> equal_range_impl(comparable_key_t const & key)
    {
        auto rng = std::ranges::equal_range(_elements, key, compare_t{});
        auto begin_pos = std::ranges::distance(_elements.begin(), rng.begin());
        auto end_pos = std::ranges::distance(_elements.begin(), rng.end());
        return {iterator{std::addressof(_elements), begin_pos},
                iterator{std::addressof(_elements), end_pos}};
    }

    template <typename comparable_key_t>
    std::pair<const_iterator, const_iterator> equal_range_impl(comparable_key_t const & key) const
    {
        auto rng = std::ranges::equal_range(_elements, key, compare_t{});
        auto begin_pos = std::ranges::distance(_elements.begin(), rng.begin());
        auto end_pos = std::ranges::distance(_elements.begin(), rng.end());
        return {const_iterator{std::addressof(_elements), begin_pos},
                const_iterator{std::addressof(_elements), end_pos}};
    }

    template <typename comparable_key_t>
    iterator lower_bound_impl(comparable_key_t const & key)
    {
        auto pos = std::ranges::distance(_elements.begin(), std::ranges::lower_bound(_elements, key, compare_t{}));
        return iterator{std::addressof(_elements), pos};
    }

    template <typename comparable_key_t>
    const_iterator lower_bound_impl(comparable_key_t const & key) const
    {
        auto pos = std::ranges::distance(_elements.begin(), std::ranges::lower_bound(_elements, key, compare_t{}));
        return const_iterator{std::addressof(_elements), pos};
    }

    template <typename comparable_key_t>
    iterator upper_bound_impl(comparable_key_t const & key)
    {
        if (_elements.size() == 1) {
            return begin();
        }
        auto pos = std::ranges::distance(_elements.begin(), std::ranges::upper_bound(_elements, key, compare_t{}));
        return iterator{std::addressof(_elements), pos};
    }

    template <typename comparable_key_t>
    const_iterator upper_bound_impl(comparable_key_t const & key) const
    {
        if (_elements.size() == 1) {
            return begin();
        }
        auto pos = std::ranges::distance(_elements.begin(), std::ranges::upper_bound(_elements, key, compare_t{}));
        return const_iterator{std::addressof(_elements), pos};
    }
};

template <std::semiregular key_t, typename compare_t>
template <bool is_const>
class sorted_vector<key_t, compare_t>::bi_iterator
{
private:
    using data_type = std::conditional_t<is_const, container_t const, container_t>;
    using data_iterator = std::ranges::iterator_t<data_type>;
    using offset_type = std::ranges::range_difference_t<data_type>;

    friend sorted_vector;

    template <bool>
    friend class bi_iterator;

    data_type * _container{};
    offset_type _index{};

    constexpr explicit bi_iterator(data_type * container, offset_type index) noexcept :
        _container{container},
        _index{index}
    {}

public:

    using value_type = std::ranges::range_value_t<data_type>;
    using reference = std::ranges::range_reference_t<data_type>;
    using difference_type = std::ranges::range_difference_t<data_type>;
    using pointer = std::conditional_t<is_const, typename container_t::const_pointer, typename container_t::pointer>;
    using iterator_category = std::random_access_iterator_tag;

    constexpr bi_iterator() = default;

    constexpr bi_iterator(bi_iterator<!is_const> iter) noexcept requires is_const :
        _container{iter._container},
        _index{iter._index}
    {}

    data_iterator base() const
    {
        return std::ranges::next(_container->begin(), _index);
    }

    constexpr reference operator*() const noexcept
    {
        assert(_container != nullptr);
        return _container->operator[](_index);
    }

    constexpr pointer operator->() const noexcept
    {
        assert(_container != nullptr);
        return std::addressof(*this);
    }

    constexpr bi_iterator & operator++() noexcept
    {
        ++_index;
        return *this;
    }

    constexpr bi_iterator operator++(int) noexcept
    {
        iterator tmp{*this};
        ++(*this);
        return tmp;
    }

    constexpr bi_iterator & operator+=(difference_type const offset) noexcept
    {
        _index += offset;
        return *this;
    }

    constexpr bi_iterator operator+(difference_type const offset) const noexcept
    {
        bi_iterator tmp{*this};
        return tmp += offset;
    }

    friend constexpr bi_iterator operator+(difference_type const offset, bi_iterator const & rhs) noexcept
    {
        return rhs + offset;
    }

    constexpr bi_iterator & operator--() noexcept
    {
        --_index;
        return *this;
    }

    constexpr bi_iterator operator--(int) noexcept
    {
        iterator tmp{*this};
        --(*this);
        return tmp;
    }

    constexpr bi_iterator & operator-=(difference_type const offset) noexcept
    {
        _index -= offset;
        return *this;
    }

    constexpr bi_iterator operator-(difference_type const offset) const noexcept
    {
        bi_iterator tmp{*this};
        return tmp -= offset;
    }

    template <bool other_const>
    constexpr difference_type operator-(bi_iterator<other_const> const & rhs) const noexcept
    {
        return _index - rhs._index;
    }

    template <bool other_const>
    constexpr bool operator==(bi_iterator<other_const> const & rhs) const noexcept
    {
        return _index == rhs._index;
    }

    template <bool other_const>
    constexpr std::strong_ordering operator<=>(bi_iterator<other_const> const & rhs) const noexcept
    {
        return _index <=> rhs._index;
    }
};

}  // namespace libjst
