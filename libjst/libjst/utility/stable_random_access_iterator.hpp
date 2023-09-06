// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a stable random access iterator tracking the index of the referenced container element.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>

namespace libjst
{

    template <std::ranges::random_access_range container_t>
    class stable_random_access_iterator {

        template <std::ranges::random_access_range>
        friend class stable_random_access_iterator;

        using base_iterator = std::ranges::iterator_t<container_t>;

    public:

        using value_type = std::iter_value_t<base_iterator>;
        using reference = std::iter_reference_t<base_iterator>;
        using difference_type = std::iter_difference_t<base_iterator>;
        using pointer = typename std::iterator_traits<base_iterator>::pointer;
        using iterator_category = typename std::iterator_traits<base_iterator>::iterator_category;

        constexpr stable_random_access_iterator() = default;
        constexpr explicit stable_random_access_iterator(container_t * base, difference_type position) noexcept :
            _base{base},
            _position{position}
        {}

        template <typename other_container_t>
            requires (std::same_as<std::remove_const_t<other_container_t>, std::remove_const_t<container_t>> &&
                      std::is_const_v<container_t> && !std::is_const_v<other_container_t>)
        constexpr stable_random_access_iterator(stable_random_access_iterator<other_container_t> other)
            noexcept :
            _base{other._base},
            _position{other._position}
        {}

        constexpr base_iterator base() const noexcept {
            return std::ranges::next(std::ranges::begin(*_base), _position);
        }

        constexpr reference operator*() const noexcept {
            return *base();
        }

        constexpr pointer operator->() const noexcept {
            return base().operator->();
        }

        constexpr reference operator[](difference_type const step) const noexcept {
            return *(base() + step);
        }

        constexpr stable_random_access_iterator & operator++() noexcept {
            ++_position;
            return *this;
        }

        constexpr stable_random_access_iterator operator++(int) noexcept {
            stable_random_access_iterator tmp{*this};
            this->operator++();
            return tmp;
        }

        constexpr stable_random_access_iterator & operator+=(difference_type const step) noexcept {
            _position += step;
            return *this;
        }

        constexpr stable_random_access_iterator & operator--() noexcept {
            --_position;
            return *this;
        }

        constexpr stable_random_access_iterator operator--(int) noexcept {
            stable_random_access_iterator tmp{*this};
            this->operator--();
            return tmp;
        }

        constexpr stable_random_access_iterator & operator-=(difference_type const step) noexcept {
            _position -= step;
            return *this;
        }

    private:

        constexpr friend stable_random_access_iterator
        operator+(stable_random_access_iterator const & lhs, difference_type const step) noexcept {
            stable_random_access_iterator tmp{lhs};
            return tmp += step;
        }

        constexpr friend stable_random_access_iterator
        operator+(difference_type const step, stable_random_access_iterator const & rhs) noexcept {
            return rhs + step;
        }

        constexpr friend stable_random_access_iterator
        operator-(stable_random_access_iterator const & lhs, difference_type const step) noexcept {
            stable_random_access_iterator tmp{lhs};
            return tmp -= step;
        }

        constexpr friend difference_type
        operator-(stable_random_access_iterator const & lhs, stable_random_access_iterator const & rhs) noexcept {
            return lhs._position - rhs._position;
        }

        constexpr friend bool
        operator==(stable_random_access_iterator const & lhs, stable_random_access_iterator const & rhs) noexcept {
            return lhs._position == rhs._position;
        }

        constexpr friend std::strong_ordering
        operator<=>(stable_random_access_iterator const & lhs, stable_random_access_iterator const & rhs) noexcept {
            return lhs._position <=> rhs._position;
        }

        container_t * _base{};
        difference_type _position{};
    };
}  // namespace libjst
