// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides common sequence interface for different sequence variant encodings.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>
#include <span>
#include <variant>

#include <seqan3/utility/detail/multi_invocable.hpp>

#include <libjst/rcms/packed_breakend_key.hpp>

namespace libjst
{
    template <std::ranges::contiguous_range source_t>
    class delta_sequence_variant : public std::ranges::view_base {
        using value_t = std::ranges::range_value_t<source_t>;
        using span_type = std::span<value_t const>;

        // class iterator_impl;
        // class sentinel_impl;

        span_type _span{};
        // iterator_impl _it{};

    public:

        using iterator = std::ranges::iterator_t<span_type>;

        constexpr delta_sequence_variant() noexcept = default;
        // constexpr delta_sequence_variant(delta_sequence_variant const & other) noexcept : _it{other._it}
        // {}
        explicit constexpr delta_sequence_variant(value_t const & snv) noexcept : _span{std::addressof(snv), 1}
        {}
        constexpr delta_sequence_variant(source_t const & insertion) noexcept : _span{insertion}
        {}

        // constexpr delta_sequence_variant & operator=(delta_sequence_variant const & other) noexcept
        // {
        //     _it = other._it;
        //     return *this;
        // }
        // ~delta_sequence_variant() noexcept = default;

        constexpr iterator begin() const noexcept {
            return _span.begin();
        }

        constexpr iterator end() const noexcept {
            return _span.end();
        }
    };

    // template <std::ranges::contiguous_range source_t>
    // class delta_sequence_variant<source_t>::iterator_impl {
    // private:
    //     friend delta_sequence_variant;
    //     friend sentinel_impl;

    //     value_t const * _it{};
    //     value_t const * _sent{};

    //     // explicit constexpr iterator_impl(value_t const & snv) noexcept : _it{_snv_table[(int)snv]}, _sent{_it + 1};
    //     // {}

    //     explicit constexpr iterator_impl(span_type indel) noexcept : _it{indel.data()}, _sent{_it + indel.size()}
    //     {}

    // public:

    //     using value_type = value_t;
    //     using reference = value_type const &;
    //     using difference_type = std::ptrdiff_t;
    //     using pointer = value_type const *;
    //     using iterator_category = std::random_access_iterator_tag;
    //     using iterator_concept = std::contiguous_iterator_tag;

    //     constexpr iterator_impl() = default;
    //     // constexpr iterator_impl(iterator_impl const & other) noexcept
    //     // {
    //     //     // when we copy the SNV, we are not comparable anymore!
    //     //     // we can only compare to itself.
    //     //     if (other._snv.has_value()) {
    //     //         _snv = other._snv;
    //     //         _sent = std::addressof(_snv.value()) + 1;
    //     //         _it = _sent - (other._sent - other._it);
    //     //     } else {
    //     //         _it = other._it;
    //     //         _sent = other._sent;
    //     //     }
    //     // }
    //     // constexpr iterator_impl & operator=(iterator_impl const & other) noexcept
    //     // {
    //     //     if (other._snv.has_value()) {
    //     //         _snv = other._snv;
    //     //         _sent = std::addressof(_snv.value()) + 1;
    //     //         _it = _sent - (other._sent - other._it);
    //     //     } else {
    //     //         _it = other._it;
    //     //         _sent = other._sent;
    //     //     }
    //     //     return *this;
    //     // }
    //     // ~iterator_impl() = default;

    //     constexpr reference operator*() const noexcept {
    //         return *_it;
    //     }

    //     constexpr pointer operator->() const noexcept {
    //         return _it;
    //     }

    //     constexpr reference operator[](difference_type const step) const noexcept {
    //         return *(_it + step);
    //     }

    //     constexpr iterator_impl & operator++() noexcept {
    //         ++_it;
    //         return *this;
    //     }

    //     constexpr iterator_impl operator++(int) noexcept {
    //         iterator_impl tmp{*this};
    //         ++(*this);
    //         return tmp;
    //     }

    //     constexpr iterator_impl & operator+=(difference_type const step) noexcept {
    //         _it += step;
    //         return *this;
    //     }

    //     constexpr iterator_impl & operator--() noexcept {
    //         --_it;
    //         return *this;
    //     }

    //     constexpr iterator_impl operator--(int) noexcept {
    //         iterator_impl tmp{*this};
    //         --(*this);
    //         return tmp;
    //     }

    //     constexpr iterator_impl & operator-=(difference_type const step) noexcept {
    //         _it -= step;
    //         return *this;
    //     }

    // private:

    //     constexpr friend iterator_impl operator+(iterator_impl lhs, difference_type const step) noexcept {
    //         return lhs += step;
    //     }

    //     constexpr friend iterator_impl operator+(difference_type const step, iterator_impl rhs) noexcept {
    //         return rhs + step;
    //     }

    //     constexpr friend iterator_impl operator-(iterator_impl lhs, difference_type const step) noexcept {
    //         return lhs -= step;
    //     }

    //     constexpr friend difference_type operator-(iterator_impl const & lhs, iterator_impl const & rhs) noexcept {
    //         // if (lhs._snv.has_value())  {
    //         //     assert(rhs._snv.has_value());
    //         //     assert(lhs._snv == rhs._snv);
    //         //     // same range different memory address
    //         //     return (rhs._sent - rhs._it) - (lhs._sent - lhs._it);
    //         // } else
    //         {
    //             return lhs._it - rhs._it;
    //         }
    //     }

    //     constexpr friend difference_type operator-(sentinel_impl const &, iterator_impl const & rhs) noexcept {
    //         return rhs._sent - rhs._it;
    //     }

    //     constexpr friend difference_type operator-(iterator_impl const & lhs, sentinel_impl const &) noexcept {
    //         return lhs._it - lhs._sent;
    //     }

    //     constexpr friend bool operator==(iterator_impl const & lhs, iterator_impl const & rhs) noexcept {
    //         return (lhs <=> rhs) == 0;
    //     }

    //     constexpr friend bool operator==(iterator_impl const & lhs, sentinel_impl const &) noexcept {
    //         return lhs._it == lhs._sent;
    //     }

    //     constexpr friend std::strong_ordering operator<=>(iterator_impl const & lhs, iterator_impl const & rhs) noexcept {
    //         // if (lhs._snv.has_value()) {
    //         //     assert(rhs._snv.has_value());
    //         //     assert(lhs._snv == rhs._snv);
    //         //     // same range different memory address
    //         //     return (lhs._sent - lhs._it) <=> (rhs._sent - rhs._it);
    //         // } else
    //         {
    //             return lhs._it <=> rhs._it;
    //         }
    //     }
    // };

    // template <std::ranges::contiguous_range source_t>
    // class delta_sequence_variant<source_t>::sentinel_impl {
    // public:
    //     sentinel_impl() = default;
    // };

}  // namespace libjst

namespace std::ranges {

    //!\cond
    template <typename source_t>
    inline constexpr bool enable_borrowed_range<libjst::delta_sequence_variant<source_t>> = true;
    //!\endcond
} // namespace std
