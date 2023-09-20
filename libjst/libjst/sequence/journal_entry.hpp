// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides internal value type of journal class.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <compare>
#include <concepts>
#include <iterator>
#include <optional>
#include <ranges>

namespace libjst
{

    /**
     * @brief Implements the value type of the journal.
     *
     * @tparam sequence_size_t The type of the position in the sequence. Must model std::integral.
     * @tparam sequence_iterator_t The type of the sequence iterator. Must model std::random_access_iterator.
     *
     * A journal entry is a key-value-pair, where the key is the sequence position and the mapped value is a
     * std::ranges::subrange over the sequence.
     * The journal entry is weakly comparable and can be ordered by the begin position of the referenced sequence in
     * the journaled sequence.
     */
    template <std::integral sequence_size_t, std::random_access_iterator sequence_iterator_t>
    class journal_entry
    {
        // ----------------------------------------------------------------------------
        // Types
        // ----------------------------------------------------------------------------
    public:
        /** @brief The type of positions in sequence space.*/
        using size_type = sequence_size_t;
        /** @brief The sequence type is implemented as std::ranges::subrange.*/
        using sequence_type = std::ranges::subrange<sequence_iterator_t, sequence_iterator_t>;

        // ----------------------------------------------------------------------------
        // Constructors, assignment and destructor
        // ----------------------------------------------------------------------------
    public:
        /**
         * @brief Construct a new journal entry object.
         *
         * The position will be zero initialized and the subrange will be default constructed.
         * This constructor is only available if the sequence type is std::default_initializable.
         */
        journal_entry()
            requires std::default_initializable<sequence_type>
        = default;

        /**
         * @brief Constructs a journal entry from a position and a sequence.
         *
         * @param position The start position of the referenced sequence segment in the journaled sequence.
         * @param segment The sequence segment referenced by this journal entry.
         */
        journal_entry(size_type position, sequence_type segment) : _position{position}, _segment{std::move(segment)}
        {
        }

        // ----------------------------------------------------------------------------
        // Accessor functions
        // ----------------------------------------------------------------------------
    public:
        /**
         * @brief Returns the begin position of the referenced segment in the journaled sequence.
         *
         * @return size_type
         * @throw Never throws.
         */
        size_type begin_position() const noexcept
        {
            return _position;
        }

        /**
         * @brief Returns the end position of the referenced segment in the journaled sequence.
         *
         * @return size_type
         * @throw Never throws.
         */
        size_type end_position() const noexcept
        {
            return begin_position() + segment().size();
        }

        /**
         * @brief Returns the referenced sequence segment.
         *
         * @return sequence_type
         * @throw Never throws if the sequence type is nothrow copy constructible.
         */
        sequence_type segment() const noexcept(std::is_nothrow_copy_constructible_v<sequence_type>)
        {
            return _segment;
        }

        // ----------------------------------------------------------------------------
        // Comparison and ordering operators
        // ----------------------------------------------------------------------------
    private:
        /**
         * @brief Equality comparison operator.
         *
         * @param lhs The left hand side of the comparison.
         * @param rhs The right hand side of the comparison.
         *
         * This operator only compares the begin positions of the journal entries to determine their equality.
         *
         * @return true if the begin positions of the journal entries are equal.
         * @return false if the begin positions of the journal entries are not equal.
         * @throws Never throws.
         */
        friend bool operator==(journal_entry const &lhs, journal_entry const &rhs) noexcept
        {
            return lhs.begin_position() == rhs.begin_position() && lhs.segment().begin() == rhs.segment().begin();
        }

        /**
         * @brief Three-way comparison operator.
         *
         * @param lhs The left hand side of the comparison.
         * @param rhs The right hand side of the comparison.
         *
         * This operator only compares the begin positions of the journal entries to determine their order.
         *
         * @returns std::weak_ordering::less if the begin position of `lhs` is less than the begin position of `rhs`.
         * @returns std::weak_ordering::greater if the begin position of `lhs` is greater than the begin position of `rhs`.
         * @returns std::weak_ordering::equivalent if the begin position of `lhs` is equal to the begin position of `rhs`.
         * @throws Never throws.
         */
        friend std::weak_ordering operator<=>(journal_entry const &lhs, journal_entry const &rhs) noexcept
        {
            return lhs.begin_position() <=> rhs.begin_position();
        }

        /**
         * @brief Tests whether the given position is covered by the journal entry.
         *
         * @param entry The journal entry whose begin and end position is used for the test.
         * @param position The position to test.
         * @return true if `position` is element of the interval `[entry.begin_position(), entry.end_position())`.
         * @return false otherwise.
         * @throws Never throws.
         */
        friend bool position_is_covered_by(journal_entry const & entry, size_type const & position) noexcept
        {
            return entry.begin_position() <= position && position < entry.end_position();
        }

        friend auto split_at(journal_entry const & entry, sequence_iterator_t split_it) noexcept
            -> std::pair<std::optional<journal_entry>, journal_entry>
        {
            if (split_it == entry.segment().begin())
                return {std::nullopt, std::move(entry)};

            size_type const split_position = entry.begin_position() +
                                                 std::ranges::distance(entry.segment().begin(), split_it);

            return {journal_entry{entry.begin_position(), sequence_type{entry.segment().begin(), split_it}},
                    journal_entry{split_position, sequence_type{split_it, entry.segment().end()}}};
        }

        // ----------------------------------------------------------------------------
        // Members
        // ----------------------------------------------------------------------------
    private:
        /** @brief The begin position of the referenced segment in the journaled sequence. */
        size_type _position{};
        /** @brief The referenced segment in the journaled sequence. */
        sequence_type _segment{};
    };

    /**
     * @brief Deduction guide for journal_entry.
     *
     * @tparam position_t The type of the position in the sequence. Must model std::integral.
     * @tparam segment_t The type of the sequence segment. Must model std::ranges::viewable_range, std::ranges::random_access_range and std::ranges::sized_range.
     */
    template <std::integral position_t, std::ranges::viewable_range segment_t>
        requires std::ranges::random_access_range<segment_t> &&
                 std::ranges::sized_range<segment_t>
    journal_entry(position_t, segment_t &&) -> journal_entry<position_t, std::ranges::iterator_t<segment_t>>;

} // namespace libjst
