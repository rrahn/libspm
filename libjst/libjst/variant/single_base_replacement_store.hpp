// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides single base replacement store.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iterator>

#include <seqan3/alphabet/concept.hpp>

#include <libcontrib/std/tag_invoke.hpp>
#include <libjst/variant/concept.hpp>

namespace libjst
{

    template <seqan3::alphabet alphabet_t>
    class single_base_replacement_store {
    private:
        template <bool>
        class iterator_type;

        class element_type;

        using data_type = std::vector<element_type>;

        data_type _data{};
    public:

        using value_type = element_type;
        using difference_type = std::ranges::range_difference_t<data_type>;
        using size_type = typename data_type::size_type;
        using iterator = iterator_type<false>;
        using const_iterator = iterator_type<true>;

        using reference = std::iter_reference_t<iterator>;
        using const_reference = std::iter_reference_t<const_iterator>;

        single_base_replacement_store() = default;

        // ----------------------------------------------------------------------------
        // Access
        // ----------------------------------------------------------------------------

        constexpr reference operator[](difference_type const offset) noexcept {
            assert(static_cast<std::size_t>(offset) < size());
            return begin()[offset];
        }

        constexpr const_reference operator[](difference_type const offset) const noexcept {
            assert(static_cast<std::size_t>(offset) < size());
            return begin()[offset];
        }

        // ----------------------------------------------------------------------------
        // Capacity
        // ----------------------------------------------------------------------------

        constexpr size_type capacity() const noexcept {
            return _data.capacity();
        }

        constexpr size_type size() const noexcept {
            return _data.size();
        }

        constexpr void reserve(size_type new_capacity) {
            _data.reserve(new_capacity);
        }

        constexpr void resize(size_type new_size) {
            _data.resize(new_size);
        }

        // ----------------------------------------------------------------------------
        // Modifier
        // ----------------------------------------------------------------------------

        constexpr void push_back(value_type value) {
            _data.push_back(std::move(value));
        }

        constexpr void emplace_back(alphabet_t value) {
            _data.emplace_back(std::move(value));
        }

        // ----------------------------------------------------------------------------
        // Iterator
        // ----------------------------------------------------------------------------

        constexpr iterator begin() noexcept { return iterator{_data.begin()}; }
        constexpr const_iterator begin() const noexcept { return const_iterator{_data.begin()}; }
        constexpr iterator end() noexcept { return iterator{_data.end()}; }
        constexpr const_iterator end() const noexcept { return const_iterator{_data.end()}; }
    };

    template <seqan3::alphabet alphabet_t>
    class single_base_replacement_store<alphabet_t>::element_type {
    private:

        alphabet_t _value{};
    public:

        element_type() = default;
        element_type(alphabet_t const & value)
            noexcept(std::is_nothrow_move_constructible_v<alphabet_t>) :
            _value{std::move(value)}
        {}

    private:

        constexpr friend size_type tag_invoke(std::tag_t<libjst::ref_span>, element_type const &) noexcept {
            return 1;
        }

        constexpr friend auto tag_invoke(std::tag_t<libjst::alt_sequence>, element_type const & me) noexcept {
            return std::views::single(me._value);
        }

        constexpr friend size_type tag_invoke(std::tag_t<libjst::effective_size>, element_type const &) noexcept {
            return 0;
        }
    };

    template <seqan3::alphabet alphabet_t>
    template <bool is_const>
    class single_base_replacement_store<alphabet_t>::iterator_type {
    private:
        template <typename t>
        using maybe_const_t = std::conditional_t<is_const, t const, t>;

        using data_iterator = std::ranges::iterator_t<maybe_const_t<data_type>>;

        data_iterator _it{};

        template <seqan3::alphabet>
        friend class single_base_replacement_store;

        template <bool>
        friend class iterator_type;

        explicit constexpr iterator_type(data_iterator it)
            noexcept(std::is_nothrow_move_constructible_v<data_iterator>) :
            _it{std::move(it)}
        {}

    public:

        using value_type = std::iter_value_t<data_iterator>;
        using reference = std::iter_reference_t<data_iterator>;
        using pointer = typename std::iterator_traits<data_iterator>::pointer;
        using difference_type = std::iter_difference_t<data_iterator>;
        using iterator_category = std::random_access_iterator_tag;
        using iterator_concept = std::contiguous_iterator_tag;

        constexpr iterator_type() = default;
        constexpr iterator_type(iterator_type<!is_const> other)
            noexcept(std::is_nothrow_move_constructible_v<data_iterator>)
            requires (is_const) :
            _it{std::move(other._it)}
        {}

        constexpr reference operator*() const noexcept {
            return *_it;
        }

        constexpr reference operator[](difference_type const offset) const noexcept {
            return *(_it + offset);
        }

        constexpr pointer operator->() const noexcept {
            return std::addressof(*_it);
        }

        constexpr iterator_type & operator++() noexcept {
            ++_it;
            return *this;
        }

        constexpr iterator_type  operator++(int) noexcept(std::is_nothrow_copy_constructible_v<data_iterator>) {
            iterator_type tmp{*this};
            ++(*this);
            return tmp;
        }

        constexpr iterator_type & operator--() noexcept {
            --_it;
            return *this;
        }

        constexpr iterator_type  operator--(int) noexcept(std::is_nothrow_copy_constructible_v<data_iterator>) {
            iterator_type tmp{*this};
            --(*this);
            return tmp;
        }

        constexpr iterator_type & operator+=(difference_type const offset) noexcept {
            _it += offset;
            return *this;
        }

        constexpr iterator_type & operator-=(difference_type const offset) noexcept {
            _it -= offset;
            return *this;
        }

    private:
        constexpr friend iterator_type operator+(iterator_type const & lhs, difference_type const offset)
            noexcept(std::is_nothrow_copy_constructible_v<data_iterator>) {
            iterator_type tmp{lhs};
            tmp += offset;
            return tmp;
        }

        constexpr friend iterator_type operator+(difference_type const offset, iterator_type const & rhs)
            noexcept(std::is_nothrow_copy_constructible_v<data_iterator>) {
            return rhs + offset;
        }

        constexpr friend iterator_type operator-(iterator_type const & lhs, difference_type const offset)
            noexcept(std::is_nothrow_copy_constructible_v<data_iterator>) {
            iterator_type tmp{lhs};
            tmp -= offset;
            return tmp;
        }

        constexpr friend difference_type operator-(iterator_type const & lhs, iterator_type const & rhs) noexcept {
            return lhs._it - rhs._it;
        }

        constexpr friend bool operator==(iterator_type const & lhs, iterator_type const & rhs) noexcept = default;
        constexpr friend auto operator<=>(iterator_type const & lhs, iterator_type const & rhs) noexcept = default;
    };

}  // namespace libjst
