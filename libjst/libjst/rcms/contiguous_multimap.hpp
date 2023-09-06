// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a map storing keys and values contiguously in two separate buffers respectively.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <bit>
#include <concepts>

#include <seqan3/core/concept/cereal.hpp>

#include <libjst/utility/sorted_vector.hpp>
#include <libjst/utility/stable_random_access_iterator.hpp>

namespace libjst
{

    template <typename key_t, typename value_t>
    struct contiguous_multimap_proxy : public std::pair<key_t, value_t>
    {
        using base_t = std::pair<key_t, value_t>;

        using base_t::base_t;

        contiguous_multimap_proxy() = default;

        template <typename _key_t, typename _value_t>
            requires std::constructible_from<key_t, _key_t &> && std::constructible_from<value_t, _value_t &>
        contiguous_multimap_proxy(contiguous_multimap_proxy<_key_t, _value_t> & other) :
            base_t{other.first, other.second}
        {}

        template <typename _key_t, typename _value_t>
            requires std::constructible_from<key_t, _key_t const &> && std::constructible_from<value_t, _value_t const &>
        contiguous_multimap_proxy(contiguous_multimap_proxy<_key_t, _value_t> const & other) :
            base_t{other.first, other.second}
        {}

        template <typename _key_t, typename _value_t>
            requires std::constructible_from<key_t, _key_t> && std::constructible_from<value_t, _value_t>
        contiguous_multimap_proxy(contiguous_multimap_proxy<_key_t, _value_t> && other) :
            base_t{std::move(other.first), std::move(other.second)}
        {}

        template <typename _key_t, typename _value_t>
            requires std::constructible_from<key_t, _key_t> && std::constructible_from<value_t, _value_t>
        contiguous_multimap_proxy(contiguous_multimap_proxy<_key_t, _value_t> const && other) :
            base_t{std::move(other.first), std::move(other.second)}
        {}

    };

    template <typename key_t, typename value_t>
    class contiguous_multimap { // contiguous_multimap
    public:

        using key_type = key_t;
        using mapped_type = value_t;

    private:

        using multiset_type = sorted_vector<key_type>;
        using data_type = std::vector<value_t>;

        template <bool is_const>
        class iterator_impl;

        multiset_type _breakends{};
        data_type _data{};

    public:

        using iterator = iterator_impl<false>;
        using const_iterator = iterator_impl<true>;
        using value_type = std::iter_value_t<iterator>;
        using size_type = typename data_type::size_type;

        constexpr contiguous_multimap() = default;

        iterator insert(value_type elem) {
            return insert_impl(std::move(elem));
        }

        iterator insert(const_iterator hint, value_type elem) {
            return insert_impl(std::move(hint), std::move(elem));
        }

        template <typename ...args_t>
            requires std::constructible_from<value_type, args_t...>
        iterator emplace(args_t &&... args) {
            return insert(value_type{(args_t &&)args...});
        }

        template <typename ...args_t>
            requires std::constructible_from<value_type, args_t...>
        iterator emplace_hint(const_iterator hint, args_t &&... args) {
            return insert(std::move(hint), value_type{(args_t &&)args...});
        }

        // void erase(const_iterator first, const_iterator last = std::ranges::next(first));
        // iterator find(key_type);

        constexpr void reserve(size_type const new_capacity) {
            _breakends.reserve(new_capacity);
            _data.reserve(new_capacity);
        }

        iterator begin() noexcept {
            return iterator{_breakends.begin(), get_data_iter(0)};
        }

        const_iterator begin() const noexcept {
            return const_iterator{_breakends.begin(), get_data_iter(0)};
        }

        iterator end() noexcept {
            return iterator{_breakends.end(), get_data_iter(std::ranges::ssize(_data))};
        }

        const_iterator end() const noexcept {
            return const_iterator{_breakends.end(), get_data_iter(std::ranges::ssize(_data))};
        }

        // ----------------------------------------------------------------------------
        // Serialisation
        // ----------------------------------------------------------------------------

        template <seqan3::cereal_input_archive archive_t>
        void load(archive_t & iarchive)
        {
            iarchive(_breakends, _data);
        }

        template <seqan3::cereal_output_archive archive_t>
        void save(archive_t & oarchive) const
        {
            oarchive(_breakends, _data);
        }

    private:

        constexpr auto get_data_iter(std::ptrdiff_t const & pos) noexcept {
            return stable_random_access_iterator{std::addressof(_data), pos};
        }

        constexpr auto get_data_iter(std::ptrdiff_t const & pos) const noexcept {
            return stable_random_access_iterator{std::addressof(_data), pos};
        }

        iterator insert_impl(value_type elem) {
            // now iterator remains stable + strong exception guarantee
            reserve(std::bit_ceil(_data.size() + 1));
            auto breakend_it = _breakends.emplace(std::move(get<0>(elem)));
            auto data_offset = std::ranges::distance(_breakends.begin(), breakend_it);
            _data.insert(std::ranges::next(_data.begin(), data_offset), std::move(get<1>(elem)));
            return iterator{std::move(breakend_it), get_data_iter(data_offset)};
        }

        iterator insert_impl(const_iterator hint, value_type elem) {
            // now iterator remains stable + strong exception guarantee
            reserve(std::bit_ceil(_data.size() + 1));
            auto [breakend_hint, data_hint] = std::move(hint).base();
            auto breakend_it = _breakends.emplace_hint(std::move(breakend_hint), std::move(get<0>(elem)));
            auto data_offset = std::ranges::distance(_breakends.begin(), breakend_it);
            _data.insert(std::ranges::next(_data.begin(), data_offset), std::move(get<1>(elem)));
            return iterator{std::move(breakend_it), get_data_iter(data_offset)};
        }
    };

    template <typename key_t, typename value_t>
    template <bool is_const>
    class contiguous_multimap<key_t, value_t>::iterator_impl {

        friend contiguous_multimap;

        template <bool>
        friend class iterator_impl;

        template <typename t>
        using maybe_const_t = std::conditional_t<is_const, t const, t>;

        using breakend_iterator = std::ranges::iterator_t<maybe_const_t<multiset_type>>;
        using data_iterator = stable_random_access_iterator<maybe_const_t<data_type>>;

        breakend_iterator _breakend_it{};
        data_iterator _data_it{};


        explicit iterator_impl(breakend_iterator breakend_it, data_iterator data_it) noexcept :
            _breakend_it{std::move(breakend_it)},
            _data_it{std::move(data_it)}
        {}

    public:

        using value_type = contiguous_multimap_proxy<std::iter_value_t<breakend_iterator> const, std::iter_value_t<data_iterator>>;
        using reference = contiguous_multimap_proxy<std::iter_value_t<breakend_iterator> const, std::iter_reference_t<data_iterator>>;
        using difference_type = std::iter_difference_t<breakend_iterator>;
        using pointer = std::optional<reference>;
        using iterator_category = std::random_access_iterator_tag;

        iterator_impl() = default;
        iterator_impl(iterator_impl<!is_const> other) noexcept requires is_const :
            _breakend_it{std::move(other._breakend_it)},
            _data_it{std::move(other._data_it)}
        {}

        constexpr std::pair<breakend_iterator, data_iterator> base() const & noexcept {
            return {_breakend_it, _data_it};
        }

        constexpr std::pair<breakend_iterator, data_iterator> base() && noexcept {
            return {std::move(_breakend_it), std::move(_data_it)};
        }

        constexpr reference operator*() const noexcept {
            return reference{*_breakend_it, *_data_it};
        }

        constexpr pointer operator->() const noexcept {
            return pointer{this->operator*()};
        }

        constexpr reference operator[](difference_type const step) const noexcept {
            return _breakend_it.operator[](step);
        }

        constexpr iterator_impl & operator++() noexcept {
            ++_breakend_it;
            ++_data_it;
            return *this;
        }

        constexpr iterator_impl operator++(int) noexcept {
            iterator_impl tmp{*this};
            ++(*this);
            return tmp;
        }

        constexpr iterator_impl & operator+=(difference_type const step) noexcept {
            _breakend_it += step;
            _data_it += step;
            return *this;
        }

        constexpr iterator_impl & operator--() noexcept {
            --_breakend_it;
            --_data_it;
            return *this;
        }

        constexpr iterator_impl operator--(int) noexcept {
            iterator_impl tmp{*this};
            --(*this);
            return tmp;
        }

        constexpr iterator_impl & operator-=(difference_type const step) noexcept {
            _breakend_it -= step;
            _data_it -= step;
            return *this;
        }

    private:

        constexpr friend iterator_impl operator+(iterator_impl lhs, difference_type const step) noexcept {
            return lhs += step;
        }

        constexpr friend iterator_impl operator+(difference_type const step, iterator_impl rhs) noexcept {
            return rhs + step;
        }

        constexpr friend iterator_impl operator-(iterator_impl lhs, difference_type const step) noexcept {
            return lhs -= step;
        }

        constexpr friend difference_type operator-(iterator_impl const & lhs, iterator_impl const & rhs) noexcept {
            return lhs._breakend_it - rhs._breakend_it;
        }

        constexpr friend bool operator==(iterator_impl const &, iterator_impl const & ) = default;
        constexpr friend auto operator<=>(iterator_impl const &, iterator_impl const & ) = default;
    };
}  // namespace libjst

namespace std {

    template <typename key1_t, typename value1_t,
              typename key2_t, typename value2_t,
              template <typename> typename T1Q,
              template <typename> typename T2Q>
        requires requires {
            typename std::common_reference_t<T1Q<key1_t>, T2Q<key2_t>>;
            typename std::common_reference_t<T1Q<value1_t>, T2Q<value2_t>>;
        }
    struct basic_common_reference<libjst::contiguous_multimap_proxy<key1_t, value1_t>,
                                  libjst::contiguous_multimap_proxy<key2_t, value2_t>,
                                  T1Q,
                                  T2Q>
    {
        using type = libjst::contiguous_multimap_proxy<
                        std::common_reference_t<T1Q<key1_t>, T2Q<key2_t>>,
                        std::common_reference_t<T1Q<value1_t>, T2Q<value2_t>>
        >;
    };

} // namespace std
