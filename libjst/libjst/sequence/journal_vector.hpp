// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a journal implementation using contiguous memory to store the elements.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <iterator>
#include <vector>
#include <ranges>
#include <type_traits>

#include <libjst/sequence/journal_entry.hpp>
#include <libjst/sequence/journal_position.hpp>

namespace libjst
{

    template <
        std::forward_iterator sequence_iterator_t,
        std::integral sequence_size_t = std::make_unsigned_t<std::iter_difference_t<sequence_iterator_t>>>
    class journal_vector
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
        using entry_type = journal_entry<sequence_size_t, sequence_iterator_t>;
        using journal_type = std::vector<entry_type>;
        using sequence_type = typename entry_type::sequence_type;

    public:
        using size_type = typename entry_type::size_type;
        using iterator = iterator_impl<false>;
        using const_iterator = iterator_impl<true>;
        using value_type = std::iter_value_t<iterator>;
        using reference = std::iter_reference_t<iterator>;
        using const_reference = std::iter_reference_t<const_iterator>;
        using difference_type = std::iter_difference_t<iterator>;

        using position_type = journal_position<iterator>;
        using const_position_type = journal_position<const_iterator>;

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

        journal_vector() : _journal{}
        {
            clear();
        }

        // ----------------------------------------------------------------------------
        // Iterators

        iterator begin() noexcept
        {
            return iterator{_journal.begin()};
        }

        const_iterator begin() const noexcept
        {
            return const_iterator{_journal.begin()};
        }

        const_iterator cbegin() const noexcept
        {
            return begin();
        }

        iterator end() noexcept
        {
            return iterator{std::prev(_journal.end())};
        }

        const_iterator end() const noexcept
        {
            return const_iterator{std::prev(_journal.end())};
        }

        const_iterator cend() const noexcept
        {
            return end();
        }

        // ----------------------------------------------------------------------------
        // Capacity

        bool empty() const noexcept
        {
            assert(!_journal.empty());
            return _journal.size() == 1;
        }

        size_type size() const noexcept
        {
            assert(!_journal.empty());
            return _journal.size() - 1;
        }

        size_type max_size() const noexcept
        {
            return std::numeric_limits<size_type>::max() - 1;
        }

        // ----------------------------------------------------------------------------
        // Modifiers

        void clear() noexcept
        {
            _journal.clear();
            _journal.emplace_back(size_type{}, sequence_type{});
        }

        position_type record_sequence_edit(position_type low, position_type high, sequence_type segment)
            requires (!std::same_as<position_type, const_position_type>)
        {
            if (low == high && std::ranges::empty(segment))
                return high;

            return record_impl(low, high, std::move(segment));
        }

        position_type record_sequence_edit(const_position_type low, const_position_type high, sequence_type segment)
        {
            return record_sequence_edit(const_cast_position(std::move(low)), const_cast_position(std::move(high)),
                                        std::move(segment));
        }

        // ----------------------------------------------------------------------------
        // Lookup

        iterator lower_bound(size_type const sequence_position) noexcept
        {
            return std::ranges::lower_bound(*this, sequence_position, std::ranges::less{}, &entry_type::begin_position);
        }

        const_iterator lower_bound(size_type const sequence_position) const noexcept
        {
            return std::ranges::lower_bound(*this, sequence_position, std::ranges::less{}, &entry_type::begin_position);
        }

        iterator upper_bound(size_type const sequence_position) noexcept
        {
            return std::ranges::upper_bound(*this, sequence_position, std::ranges::less{}, &entry_type::begin_position);
        }

        const_iterator upper_bound(size_type const sequence_position) const noexcept
        {
            return std::ranges::upper_bound(*this, sequence_position, std::ranges::less{}, &entry_type::begin_position);
        }

        position_type find(size_type const sequence_position) noexcept
        {
            if (sequence_position >= _journal.back().end_position())
                return position_type{end()};

            auto journal_it = upper_bound(sequence_position) - 1;
            assert(journal_it->position_is_covered_by(sequence_position));
            difference_type const segment_offset = sequence_position - journal_it->begin_position();
            auto segment_it = std::ranges::next(journal_it->segment().begin(), segment_offset);
            return position_type{std::move(journal_it), std::move(segment_it)};
        }

        const_position_type find(size_type const sequence_position) const noexcept
        {
            if (sequence_position >= _journal.back().end_position())
                return const_position_type{end()};

            auto journal_it = upper_bound(sequence_position) - 1;
            assert(journal_it->position_is_covered_by(sequence_position));
            difference_type const segment_offset = sequence_position - journal_it->begin_position();
            auto segment_it = std::ranges::next(journal_it->segment().begin(), segment_offset);
            return const_position_type{std::move(journal_it), std::move(segment_it)};
        }

        // ----------------------------------------------------------------------------
        // Utility functions
    private:

        position_type const_cast_position(const_position_type const & pos) noexcept
        {
            auto [const_journal_it, sequence_it] = pos;
            return position_type{begin() + (const_journal_it - cbegin()), sequence_it};
        }

        position_type record_impl(position_type low, position_type high, sequence_type new_segment)
        {
            assert(low <= high);

            difference_type const deletion_size = to_sequence_position(high) - to_sequence_position(low);
            difference_type const insertion_size = std::ranges::ssize(new_segment);

            std::array<entry_type, 2> entries_marked_for_insertion{};
            size_t marked_insert_entries{};

            auto [maybe_low_prefix, low_suffix] = split_at(low);
            auto [ignore_high_prefix, high_suffix] = split_at(high);

            if (!std::ranges::empty(maybe_low_prefix.segment()))
                entries_marked_for_insertion[marked_insert_entries++] = std::move(maybe_low_prefix);

            if (insertion_size > 0)
                entries_marked_for_insertion[marked_insert_entries++] =
                        entry_type{low_suffix.begin_position(), std::move(new_segment)};

            *(high.journal_it) = std::move(high_suffix);
            auto insert_position = _journal.erase(low.journal_it.base(), high.journal_it.base());
            auto insert_entries_it = std::make_move_iterator(entries_marked_for_insertion.begin());
            auto first_inserted = _journal.insert(insert_position,
                                                  insert_entries_it,
                                                  insert_entries_it + marked_insert_entries);

            update_positions_of_remaining_entries(first_inserted + marked_insert_entries, insertion_size - deletion_size);

            assert(check_journal_invariants());

            return position_type{iterator{first_inserted + marked_insert_entries - (insertion_size > 0)}};
        }

        void update_positions_of_remaining_entries(journal_type::iterator journal_it, difference_type offset) noexcept
        {
            for (; journal_it != std::ranges::end(_journal); ++journal_it)
                *journal_it = entry_type{journal_it->begin_position() + offset, std::move(journal_it->segment())};
        }

        bool check_journal_invariants() const noexcept {
            // 1. check if first entry starts at position 0
            // 2. check if all adjacent entries are non-overlapping such that end position of smaller entry is equal to begin position of successor
            if (_journal.empty() || _journal.begin()->begin_position() != 0)
                return false;

            for (auto it = begin(); it != end(); ++it)
            {
                auto const &entry = *it;
                auto const &next_entry = *std::ranges::next(it);
                if (entry.end_position() != next_entry.begin_position())
                    return false;
            }
            return true;
        }
    };

    template <std::forward_iterator sequence_iterator_t, std::integral sequence_size_t>
    template <bool is_const>
    class journal_vector<sequence_iterator_t, sequence_size_t>::iterator_impl
    {
        friend journal_vector;

        template <bool>
        friend class iterator_impl;

        // ----------------------------------------------------------------------------
        // Member Types
        // ----------------------------------------------------------------------------
    private:
        using maybe_const_journal_type = std::conditional_t<is_const, journal_type const, journal_type>;
        using journal_iterator = std::ranges::iterator_t<maybe_const_journal_type>;

    public:
        using value_type = std::iter_value_t<journal_iterator>;
        using reference = std::iter_reference_t<journal_iterator>;
        using difference_type = std::iter_difference_t<journal_iterator>;
        using pointer = typename std::iterator_traits<journal_iterator>::pointer;
        using iterator_category = std::random_access_iterator_tag;
        using iterator_concept = std::contiguous_iterator_tag;

        // ----------------------------------------------------------------------------
        // Member Variables
        // ----------------------------------------------------------------------------
    private:

        journal_iterator _journal_it{};

        // ----------------------------------------------------------------------------
        // Member Functions
        // ----------------------------------------------------------------------------

        // ----------------------------------------------------------------------------
        // Constructors, assignment and destructor
    private:

        explicit iterator_impl(journal_iterator journal_it) :
            _journal_it{std::move(journal_it)}
        {
        }

    public:

        iterator_impl() = default;

        iterator_impl(iterator_impl<!is_const> other) requires
            (is_const && std::convertible_to<std::ranges::iterator_t<journal_type>, journal_iterator>) :
            _journal_it{std::move(other._journal_it)}
        {
        }

        // ----------------------------------------------------------------------------
        // Element access

        reference operator*() const noexcept
        {
            return *_journal_it;
        }

        pointer operator->() const noexcept
        {
            return std::addressof(this->operator*());
        }

        reference operator[](difference_type const offset) const noexcept
        {
            return *(*this + offset);
        }

        journal_iterator base() const & noexcept
        {
            return _journal_it;
        }

        journal_iterator base() && noexcept
        {
            return std::move(_journal_it);
        }

        // ----------------------------------------------------------------------------
        // Arithmetic operators

        iterator_impl &operator++() noexcept
        {
            ++_journal_it;
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
            _journal_it += count;
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
            --_journal_it;
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
            _journal_it -= count;
            return *this;
        }

        friend iterator_impl operator-(iterator_impl const &lhs, difference_type const rhs) noexcept
        {
            iterator_impl tmp{lhs};
            tmp -= rhs;
            return tmp;
        }

        friend difference_type operator-(iterator_impl const &lhs, iterator_impl const &rhs) noexcept
        {
            return lhs._journal_it - rhs._journal_it;
        }

        // ----------------------------------------------------------------------------
        // Comparison operators

        friend bool operator==(iterator_impl const &lhs, iterator_impl const &rhs) noexcept = default;

        friend std::strong_ordering operator<=>(iterator_impl const &lhs, iterator_impl const &rhs) noexcept = default;
    };

}  // namespace libjst
