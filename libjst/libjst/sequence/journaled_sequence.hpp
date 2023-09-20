// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of libjst::journaled_sequence.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <iterator>
#include <limits>
#include <ranges>

#include <libjst/sequence/journal_vector.hpp>

namespace libjst
{

    template <
        std::ranges::random_access_range sequence_t,
        std::integral sequence_size_t = std::ranges::range_size_t<sequence_t>>
    class journaled_sequence
    {
        // ----------------------------------------------------------------------------
        // Member classes
        // ----------------------------------------------------------------------------
    public:
        template <bool>
        class iterator_impl;

        // ----------------------------------------------------------------------------
        // Member Types
        // ----------------------------------------------------------------------------
    private:
        using journal_type = journal_vector<std::ranges::iterator_t<sequence_t>, sequence_size_t>;
        using entry_type = std::ranges::range_value_t<journal_type>;

    public:
        using size_type = typename journal_type::size_type;
        using iterator = iterator_impl<false>;
        using const_iterator = iterator_impl<true>;
        using value_type = std::iter_value_t<iterator>;
        using reference = std::iter_reference_t<iterator>;
        using const_reference = std::iter_reference_t<const_iterator>;
        using difference_type = std::iter_difference_t<iterator>;
        using segment_type = typename entry_type::sequence_type;

        // ----------------------------------------------------------------------------
        // Member Variables
        // ----------------------------------------------------------------------------
    private:
        journal_type _journal{};

        // ----------------------------------------------------------------------------
        // Member functions
        // ----------------------------------------------------------------------------
    public:
        // ----------------------------------------------------------------------------
        // Constructors, assignment and destructor

        journaled_sequence() = default;

        journaled_sequence(segment_type initial_sequence)
        {
            insert(begin(), std::move(initial_sequence));
        }

        // ----------------------------------------------------------------------------
        // Iterators

        iterator begin() noexcept
        {
            return iterator{std::addressof(_journal), _journal.begin()};
        }

        const_iterator begin() const noexcept
        {
            return const_iterator{std::addressof(_journal), _journal.begin()};
        }

        iterator end() noexcept
        {
            return iterator{std::addressof(_journal), _journal.end()};
        }

        const_iterator end() const noexcept
        {
            return const_iterator{std::addressof(_journal), _journal.end()};
        }

        // ----------------------------------------------------------------------------
        // Capacity

        bool empty() const noexcept
        {
            return _journal.empty();
        }

        size_type size() const noexcept
        {
            return _journal.end()->begin_position();
        }

        size_type max_size() const noexcept
        {
            return std::numeric_limits<size_type>::max();
        }

        // ----------------------------------------------------------------------------
        // Modifiers

        void clear() noexcept
        {
            _journal.clear();
        }

        iterator insert(const_iterator position, segment_type segment)
        {
            return replace(position, position, std::move(segment));
        }

        iterator erase(const_iterator low)
        {
            return erase(low, low + 1);
        }

        iterator erase(const_iterator low, const_iterator high)
        {
            assert(low <= high);
            return replace(std::move(low), std::move(high), segment_type{});
        }

        iterator replace(const_iterator low, const_iterator high, segment_type segment)
        {
            return iterator{std::addressof(_journal),
                            _journal.record_sequence_edit(std::move(low), std::move(high), std::move(segment))};
        }
    };

    // deduction guide
    template <std::ranges::viewable_range sequence_range_t>
    journaled_sequence(sequence_range_t &&) -> journaled_sequence<sequence_range_t>;

    template <std::ranges::forward_range sequence_t, std::integral sequence_size_t>
    template <bool is_const>
    class journaled_sequence<sequence_t, sequence_size_t>::iterator_impl
    {
        friend journaled_sequence;

        template <bool>
        friend class iterator_impl;

        // ----------------------------------------------------------------------------
        // Member Types
        // ----------------------------------------------------------------------------
    private:

        using maybe_const_journal_type = std::conditional_t<is_const, journal_type const, journal_type>;
        using journal_iterator = std::ranges::iterator_t<maybe_const_journal_type>;

        using position_type = journal_position<journal_iterator>;
        using sequence_iterator = typename position_type::sequence_iterator;

    public:
        using value_type = std::iter_value_t<sequence_iterator>;
        using reference = std::iter_reference_t<sequence_iterator>;
        using difference_type = std::iter_difference_t<journal_iterator>;
        using iterator_category = std::random_access_iterator_tag;
        using pointer = typename std::iterator_traits<sequence_iterator>::pointer;

        // ----------------------------------------------------------------------------
        // Member Variables
        // ----------------------------------------------------------------------------
    private:
        maybe_const_journal_type *_journal{};
        journal_iterator _journal_it{};
        sequence_iterator _sequence_it{};

        // ----------------------------------------------------------------------------
        // Member Functions
        // ----------------------------------------------------------------------------

        // ----------------------------------------------------------------------------
        // Constructors, assignment and destructor
    private:

        explicit iterator_impl(maybe_const_journal_type *journal, position_type journal_position) :
            _journal{journal},
            _journal_it{std::move(journal_position.journal_it)},
            _sequence_it{std::move(journal_position.sequence_it)}
        {
        }

    public:

        iterator_impl() = default;

        iterator_impl(iterator_impl<!is_const> other) requires is_const :
            _journal{std::move(other._journal)},
            _journal_it{std::move(other._journal_it)},
            _sequence_it{std::move(other._sequence_it)}
        {
        }

        // ----------------------------------------------------------------------------
        // Element access

        reference operator*() const noexcept
        {
            return *_sequence_it;
        }

        pointer operator->() const noexcept
        {
            return std::addressof(this->operator*());
        }

        reference operator[](difference_type const offset) const noexcept
        {
            return *(*this + offset);
        }
    private:

        operator position_type() const & noexcept
        {
            return position_type{_journal_it, _sequence_it};
        }

        operator position_type() && noexcept
        {
            return position_type{std::move(_journal_it), std::move(_sequence_it)};
        }

        // ----------------------------------------------------------------------------
        // Arithmetic operators
    public:

        iterator_impl &operator++() noexcept
        {
            assert(_sequence_it < std::ranges::end(_journal_it->segment()));
            if (++_sequence_it == std::ranges::end(_journal_it->segment()))
            {
                ++_journal_it;
                _sequence_it = std::ranges::begin(_journal_it->segment());
            }
            return *this;
        }

        iterator_impl operator++(int) noexcept
        {
            iterator_impl tmp{*this};
            ++(*this);
            return tmp;
        }

        iterator_impl &operator+=(difference_type count) noexcept
        {
            if (auto next_position = current_position() + count; position_is_covered_by(*_journal_it, next_position))
                _sequence_it += count;
            else
                *this = iterator_impl{_journal, _journal->find(next_position)};

            return *this;
        }

        friend iterator_impl operator+(iterator_impl const &lhs, difference_type const rhs) noexcept
        {
            iterator_impl tmp{lhs};
            tmp += rhs;
            return tmp;
        }

        friend iterator_impl operator+(difference_type const lhs, iterator_impl const &rhs) noexcept
        {
            return rhs + lhs;
        }

        iterator_impl &operator--() noexcept
        {
            if (_sequence_it == std::ranges::begin(_journal_it->segment()))
            {
                --_journal_it;
                _sequence_it = std::ranges::end(_journal_it->segment());
            }
            --_sequence_it;
            return *this;
        }

        iterator_impl operator--(int) noexcept
        {
            iterator_impl tmp{*this};
            --(*this);
            return tmp;
        }

        iterator_impl &operator-=(difference_type const count) noexcept
        {
            return (*this) += -count;
        }

        friend iterator_impl operator-(iterator_impl const &lhs, difference_type const rhs) noexcept
        {
            iterator_impl tmp{lhs};
            tmp -= rhs;
            return tmp;
        }

        friend difference_type operator-(iterator_impl const &lhs, iterator_impl const &rhs) noexcept
        {
            return lhs.current_position() - rhs.current_position();
        }

        // ----------------------------------------------------------------------------
        // Comparison operators

        friend bool operator==(iterator_impl const &lhs, iterator_impl const &rhs) noexcept
        {
            return lhs._journal_it == rhs._journal_it && lhs._sequence_it == rhs._sequence_it;
        }

        friend std::strong_ordering operator<=>(iterator_impl const &lhs, iterator_impl const &rhs) noexcept
        {
            return lhs.current_position() <=> rhs.current_position();
        }

        // ----------------------------------------------------------------------------
        // Utility functions
    private:
        difference_type current_position() const noexcept
        {
            assert(static_cast<difference_type>(_journal_it->begin_position()) >= 0);

            return _journal_it->begin_position() +
                   std::ranges::distance(std::ranges::begin(_journal_it->segment()), _sequence_it);
        }
    };

} // namespace libjst
