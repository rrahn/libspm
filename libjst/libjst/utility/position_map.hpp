// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::position_map.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <compare>
#include <concepts>
#include <memory_resource>
#include <vector>

namespace libjst
{
template <std::integral key_t,
          typename mapped_t,
          typename compare_t = std::less<key_t>>
class position_map
{
    using keys_t = std::vector<key_t>;
    using mapped_values_t = std::vector<mapped_t>;
public:
    keys_t _keys; //!< The container holding the elements.
    mapped_values_t _values;

    template <bool is_const>
    class map_proxy;

    template <bool is_const>
    class map_iterator;

public:

    using value_type = std::pair<key_t const, mapped_t>;
    using iterator = map_iterator<false>;
    using const_iterator = map_iterator<true>;
    using size_type = size_t;
    using reference = map_proxy<false>;
    using const_reference = map_proxy<true>;

    position_map() noexcept = default;
    position_map(position_map const &) noexcept = default;
    position_map(position_map &&) noexcept = default;
    position_map & operator=(position_map const &) noexcept = default;
    position_map & operator=(position_map &&) noexcept = default;
    ~position_map() = default;

    /*!\name Iterators
     * \{
     */
    constexpr iterator begin() noexcept
    {
        return iterator{this, 0};
    }

    constexpr const_iterator begin() const noexcept
    {
        return const_iterator{this, 0};
    }

    constexpr iterator end() noexcept
    {
        return iterator{this, size()};
    }

    constexpr const_iterator end() const noexcept
    {
        return const_iterator{this, size()};
    }
    //!\}

    /*!\name Cpacity
     * \{
     */
    constexpr bool empty() const noexcept
    {
        return _keys.empty();
    }

    constexpr size_type size() const noexcept
    {
        return _keys.size();
    }

    constexpr size_type max_size() const noexcept
    {
        return _keys.max_size();
    }
    //!\}

    /*!\name Modifiers
     * \{
     */
    constexpr auto clear() noexcept
    {
        _keys.clear();
        _values.clear();
    }

    auto insert(value_type value) noexcept// if no_throw_move_assignable
    {
        // std::cout << "Insert at " << value.first << "\n";
        // find insert position

        if (auto it = std::ranges::lower_bound(_keys, value.first); it == std::ranges::end(_keys) || *it != value.first) {
            auto insert_position = std::ranges::distance(std::ranges::begin(_keys), it);
            // std::cout << "insert position = " << insert_position << "\n";
            _keys.emplace(std::ranges::next(std::ranges::begin(_keys), insert_position), value.first);
            _values.emplace(std::ranges::next(std::ranges::begin(_values), insert_position), std::move(value.second));
            return std::pair{iterator{this, insert_position}, true};
        }
        return std::pair{end(), false};
    }

    // iterator insert(const_iterator hint, value_type const & value)
    // {
    //     return insert_impl(hint, value);
    // }

    // iterator insert(const_iterator hint, value_type && value)
    // {
    //     return insert_impl(hint, std::move(value));
    // }

    // template <typename ...args_t>
    // iterator emplace(args_t && ...args)
    // {
    //     return insert(value_type(std::forward<args_t>(args)...));
    // }

    // template <typename ...args_t>
    // iterator emplace_hint(const_iterator hint, args_t && ...args)
    // {
    //     return insert(hint, value_type(std::forward<args_t>(args)...));
    // }

    // iterator erase(iterator pos)
    // {
    //     return iterator{_elements.erase(pos.base())};
    // }

    // iterator erase(const_iterator pos)
    // {
    //     return iterator{_elements.erase(pos.base())};
    // }

    // iterator erase(const_iterator first, const_iterator last)
    // {
    //     return iterator{_elements.erase(first.base(), last.base())};
    // }

    // size_type erase(key_t const & key)
    // {
    //     auto [first, last] = equal_range_impl(key);
    //     size_t erased_elements = std::ranges::distance(first, last);
    //     _elements.erase(first.base(), last.base());
    //     return erased_elements;
    // }
    // constexpr void swap();
    // constexpr void extract();
    // constexpr void merge();
    //!\}

    /*!\name Lookup
     * \{
     */
    // size_type count(key_t const & key) const
    // {
    //     auto [first, last] = equal_range(key);
    //     return std::ranges::distance(first, last);
    // }

    // template<typename other_key_t>
    // size_type count(other_key_t const & key) const
    // {
    //     auto [first, last] = equal_range(key);
    //     return std::ranges::distance(first, last);
    // }

    // iterator find(key_t const & key)
    // {
    //     return find_impl(key);
    // }

    // const_iterator find(key_t const & key) const
    // {
    //     return find_impl(key);
    // }

    // template<typename comparable_key_t>
    // iterator find(comparable_key_t const & key)
    // {
    //     return find_impl(key);
    // }

    // template<typename comparable_key_t>
    // const_iterator find(comparable_key_t const & key) const
    // {
    //     return find_impl(key);
    // }

    // bool contains(key_t const & key) const
    // {
    //     return find(key) != end();
    // }

    // template<typename comparable_key_t>
    // bool contains(comparable_key_t const & key) const
    // {
    //     return find(key) != end();
    // }

    // std::pair<iterator, iterator> equal_range(key_t const & key)
    // {
    //     return equal_range_impl(key);
    // }
    // std::pair<const_iterator, const_iterator> equal_range(key_t const & key) const
    // {
    //     return equal_range_impl(key);
    // }

    // template<typename other_key_t>
    // std::pair<iterator, iterator> equal_range(other_key_t const & key)
    // {
    //     return equal_range_impl(key);
    // }

    // template<typename other_key_t>
    // std::pair<const_iterator, const_iterator> equal_range(other_key_t const & key) const
    // {
    //     return equal_range_impl(key);
    // }

    // iterator lower_bound(key_t const & key)
    // {
    //     return lower_bound_impl(key);
    // }

    // const_iterator lower_bound(key_t const & key) const
    // {
    //     return lower_bound_impl(key);
    // }

    // template<typename other_key_t>
    //     requires std::strict_weak_order<compare_t, key_t, other_key_t>
    // iterator lower_bound(other_key_t const & key)
    // {
    //     return lower_bound_impl(key);
    // }

    // template<typename other_key_t>
    //     requires std::strict_weak_order<compare_t, key_t, other_key_t>
    // const_iterator lower_bound(other_key_t const & key) const
    // {
    //     return lower_bound_impl(key);
    // }

    iterator lower_bound(key_t key)
    {
        return lower_bound_impl(*this, key);
    }
    const_iterator lower_bound(key_t key) const
    {
        return lower_bound_impl(*this, key);
    }

    // template<typename other_key_t>
    // iterator upper_bound(other_key_t const & key)
    // {
    //     return upper_bound_impl(key);
    // }

    // template<typename other_key_t>
    // const_iterator upper_bound(other_key_t const & key) const
    // {
    //     return upper_bound_impl(key);
    // }

    iterator upper_bound(key_t const & key)
    {
        return upper_bound_impl(*this, key);
    }

    const_iterator upper_bound(key_t const & key) const
    {
        return upper_bound_impl(*this, key);
    }
    //!\}

    /*!\name Comparison
     * \{
     */
    bool operator==(position_map const & other) const = default;
    std::strong_ordering operator<=>(position_map const & other) const = default;
    //!\}

private:

    // template <typename comparable_key_t>
    // iterator find_impl(comparable_key_t && key)
    // {
    //     auto base_it = std::ranges::lower_bound(_elements, key, compare_t{});
    //     if (compare_t{}(key, *base_it)) // not identical
    //         return iterator{_elements.end()};

    //     return iterator{base_it};
    // }

    // template <typename comparable_key_t>
    // const_iterator find_impl(comparable_key_t && key) const
    // {
    //     auto base_it = std::ranges::lower_bound(_elements, key, compare_t{});
    //     if (compare_t{}(key, *base_it)) // not identical
    //         return const_iterator{_elements.end()};

    //     return const_iterator{base_it};
    // }

    // template <typename value_t>
    // iterator insert_impl(const_iterator hint, value_t && value)
    // {
    //     // We want to insert before the hint but only if it is ok to insert here.
    //     if (_elements.empty())
    //     {
    //         _elements.push_back(std::forward<value_t>(value));
    //         return iterator{_elements.begin()};
    //     }

    //     compare_t compare{};

    //     auto hint_base =  hint.base();
    //     bool const at_end = hint_base == _elements.end();
    //     // (at_end or value < value at hint) && value >= value directly before hint
    //     hint_base = ((at_end || compare(value, *hint_base)) && (hint_base == _elements.begin() || !compare(value, *(std::ranges::prev(hint_base)))))
    //                 ? hint_base
    //                 : std::ranges::upper_bound(_elements, value, compare);

    //     return iterator{_elements.insert(hint_base, std::forward<value_t>(value))};
    // }

    // template <typename comparable_key_t>
    // std::pair<iterator, iterator> equal_range_impl(comparable_key_t const & key)
    // {
    //     auto rng = std::ranges::equal_range(_elements, key, compare_t{});
    //     return {iterator{rng.begin()}, iterator{rng.end()}};
    // }

    // template <typename comparable_key_t>
    // std::pair<const_iterator, const_iterator> equal_range_impl(comparable_key_t const & key) const
    // {
    //     auto rng = std::ranges::equal_range(_elements, key, compare_t{});
    //     return {const_iterator{rng.begin()}, const_iterator{rng.end()}};
    // }

    template <typename map_t>
    static auto lower_bound_impl(map_t & me, key_t key) {
        using iterator_t = std::conditional_t<std::is_const_v<map_t>, const_iterator, iterator>;

        return iterator_t{&me, std::ranges::distance(std::ranges::begin(me._keys),
                                                     std::ranges::lower_bound(me._keys, key))};
    }

    // template <typename comparable_key_t>
    // const_iterator lower_bound_impl(comparable_key_t const & key) const
    // {
    //     return const_iterator{std::ranges::lower_bound(_elements, key, compare_t{})};
    // }

    template <typename map_t>
    static auto upper_bound_impl(map_t & me, key_t key) {
        using iterator_t = std::conditional_t<std::is_const_v<map_t>, const_iterator, iterator>;

        return iterator_t{&me, std::ranges::distance(std::ranges::begin(me._keys),
                                                     std::ranges::upper_bound(me._keys, key))};
    }

    // template <typename comparable_key_t>
    // iterator upper_bound_impl(comparable_key_t const & key)
    // {
    //     return (_elements.size() == 1) ? begin() : iterator{std::ranges::upper_bound(_elements, key, compare_t{})};
    // }

    // template <typename comparable_key_t>
    // const_iterator upper_bound_impl(comparable_key_t const & key) const
    // {
    //     return (_elements.size() == 1) ? begin() : const_iterator{std::ranges::upper_bound(_elements, key, compare_t{})};
    // }
};

template <std::integral key_t, typename mapped_t, typename compare_t>
template <bool is_const>
class position_map<key_t, mapped_t, compare_t>::map_proxy {

    using maybe_const_value_t = std::conditional_t<is_const, mapped_t const, mapped_t>;

public:
    using first_type = key_t;
    using second_type = maybe_const_value_t;

    first_type first{};
    second_type second{};

    map_proxy() = default;
    map_proxy(key_t key, maybe_const_value_t value) noexcept : first{key}, second{std::move(value)}
    {}

    constexpr operator value_type() noexcept {
        return value_type{first, second};
    }
};

template <std::integral key_t, typename mapped_t, typename compare_t>
template <bool is_const>
class position_map<key_t, mapped_t, compare_t>::map_iterator
{
private:
    using maybe_const_map_t = std::conditional_t<is_const, position_map const, position_map>;
    using size_type = typename position_map::size_type;

    template <bool>
    friend class map_iterator;

    maybe_const_map_t * _source{};
    size_type _position{};
public:

    using value_type = typename position_map::value_type;
    using reference = map_proxy<is_const>;
    using difference_type = std::ptrdiff_t;
    using pointer = void;
    using iterator_category = std::random_access_iterator_tag;

    constexpr map_iterator() = default;
    constexpr explicit map_iterator(maybe_const_map_t * source, size_type position) noexcept :
        _source{source},
        _position{std::move(position)}
    {}

    constexpr map_iterator(map_iterator<!is_const> iter) noexcept
        requires is_const :
        _source{iter._source},
        _position{std::move(iter._position)}
    {}

    constexpr reference operator*() const noexcept
    {
        return {_source->_keys[_position], _source->_values[_position]};
    }

    constexpr map_iterator & operator++() noexcept
    {
        ++_position;
        return *this;
    }

    constexpr map_iterator operator++(int) noexcept
    {
        iterator tmp{*this};
        ++(*this);
        return tmp;
    }

    constexpr map_iterator & operator+=(difference_type const offset) noexcept
    {
        _position += offset;
        return *this;
    }

    constexpr map_iterator operator+(difference_type const offset) const noexcept
    {
        map_iterator tmp{*this};
        return tmp += offset;
    }

    friend constexpr map_iterator operator+(difference_type const offset, map_iterator const & rhs) noexcept
    {
        return rhs + offset;
    }

    constexpr map_iterator & operator--() noexcept
    {
        --_position;
        return *this;
    }

    constexpr map_iterator operator--(int) noexcept
    {
        iterator tmp{*this};
        --(*this);
        return tmp;
    }

    constexpr map_iterator & operator-=(difference_type const offset) noexcept
    {
        _position -= offset;
        return *this;
    }

    constexpr map_iterator operator-(difference_type const offset) const noexcept
    {
        map_iterator tmp{*this};
        return tmp -= offset;
    }

    template <bool other_const>
    constexpr difference_type operator-(map_iterator<other_const> const & rhs) const noexcept
    {
        return _position - rhs._position;
    }

    template <bool other_const>
    constexpr bool operator==(map_iterator<other_const> const & rhs) const noexcept
    {
        return _position == rhs._position;
    }

    template <bool other_const>
    constexpr std::strong_ordering operator<=>(map_iterator<other_const> const & rhs) const noexcept
    {
        return _position <=> rhs._position;
    }
};

}  // namespace libjst

namespace std {

// template <integral key_t, typename mapped_t, typename compare_t, bool is_const>
// struct tuple_size<libjst::position_map<key_t, mapped_t, compare_t>::template map_proxy<is_const>>
//     : public integral_constant<2>
// {};

// template <integral key_t, typename mapped_t, typename compare_t, bool is_const>
// struct tuple_element<0, typename libjst::position_map<key_t, mapped_t, compare_t>::map_proxy<is_const>>
//     : public type_identity<typename libjst::position_map<key_t, mapped_t, compare_t>::map_proxy<is_const>::first_type>
// {};

// template <integral key_t, typename mapped_t, typename compare_t, bool is_const>
// struct tuple_element<1, typename libjst::position_map<key_t, mapped_t, compare_t>::map_proxy<is_const>>
//     : public type_identity<typename libjst::position_map<key_t, mapped_t, compare_t>::map_proxy<is_const>::second_type>
// {};
} // namespace std
