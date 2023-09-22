// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/**
 * @file
 * @brief Provides internal value type of journal class.
 * @author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
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
     * @tparam sequence_size_t The type of the position in the sequence. Must model:
     * @rst
     * |   ``std::integral``
     * @endrst
     * @tparam sequence_iterator_t The type of the sequence iterator. Must model:
     * @rst
     * |   ``std::random_access_iterator``
     * @endrst
     *
     * A journal entry is a key-value-pair, where the key is the sequence position and the mapped value is a
     * std::ranges::subrange over the sequence.
     * The journal entry is weakly comparable and can be ordered by the begin position of the referenced sequence in
     * the journaled sequence.
     */
    template <std::integral sequence_size_t, std::random_access_iterator sequence_iterator_t>
    class journal_entry
    {
        /** @name Public member types */
        ///@{
    public:
        /** @brief The type of positions in sequence space.*/
        using size_type = sequence_size_t;
        /** @brief The sequence type is implemented as std::ranges::subrange.*/
        using sequence_type = std::ranges::subrange<sequence_iterator_t, sequence_iterator_t>;
        ///@}

        /** @name Constructors, assignment and destructor */
        ///@{
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
        ///@}

        /** @name Accessors */
        ///@{
    public:
        /**
         * @brief Returns the begin position of the referenced segment in the journaled sequence.
         *
         * @returns size_type
         */
        size_type begin_position() const noexcept
        {
            return _position;
        }

        /**
         * @brief Returns the end position of the referenced segment in the journaled sequence.
         *
         * @returns size_type
         */
        size_type end_position() const noexcept
        {
            return begin_position() + segment().size();
        }

        /**
         * @brief Returns the referenced sequence segment.
         *
         * @returns sequence_type
         * @exception implementation_defined If the sequence type is nothrow copy constructible this function never throws.
         */
        sequence_type segment() const noexcept(std::is_nothrow_copy_constructible_v<sequence_type>)
        {
            return _segment;
        }
        ///@}

        /** @name Utilities */
        ///@{
        /**
         * @brief Tests whether the given position is covered by the journal entry.
         *
         * @param position The position to test.
         * @return `true` if `position` is element of the interval `[entry.begin_position(), entry.end_position())`,
         *         `false` otherwise.
         */
        bool position_is_covered_by(size_type const & position) const noexcept
        {
            return begin_position() <= position && position < end_position();
        }

        /**
         * @brief Splits the journal entry at the given position.
         *
         * @param split_it The iterator pointing to the split position.
         * @retval std::pair<journal_entry, journal_entry> A pair of journal entries, where the first entry covers the
         *         interval `[segment().begin(), split_it)` and the second entry covers the interval
         *         `[split_it, segment().end())`.
         * @rst
         * .. note:: If ``split_it`` is out of bounds, either the first entry or the second entry will have an
         *           empty segment.
         * @endrst
         *
         * @rst
         * .. danger:: If ``split_it`` is not a valid iterator of the sequence segment, the behaviour is undefined.
         * @endrst
         */
        auto split_at(sequence_iterator_t split_it) const noexcept -> std::pair<journal_entry, journal_entry>
        {
            auto split_offset = std::clamp<int32_t>(std::ranges::distance(segment().begin(), std::move(split_it)),
                                                    0, std::ranges::ssize(segment()));
            split_it = std::ranges::next(segment().begin(), split_offset);
            auto left_entry = journal_entry{begin_position(),
                                            sequence_type{segment().begin(), split_it}};
            auto right_entry = journal_entry{begin_position() + split_offset,
                                             sequence_type{split_it, segment().end()}};

            return {std::move(left_entry), std::move(right_entry)};
        }
        ///@}

        /** @name Non-member functions */
        ///@{
    public:
        /**
         * @brief Equality comparison.
         *
         * @param lhs The left hand side of the comparison.
         * @param rhs The right hand side of the comparison.
         *
         * @result `true` if both entries refer to the same segment, `false` otherwise.
         * @rst
         * .. note::
         *    This operator does not compare the segments lexicographically, but only if they point to the same
         *    memory range in the original sequence.
         * @endrst
         */
        friend bool operator==(journal_entry const &lhs, journal_entry const &rhs) noexcept
        {
            return lhs.begin_position() == rhs.begin_position() &&
                   lhs.segment().begin() == rhs.segment().begin() &&
                   lhs.end_position() == rhs.end_position();
        }

        /**
         * @brief Three-way comparison.
         *
         * @param lhs The left hand side of the comparison.
         * @param rhs The right hand side of the comparison.
         *
         * This operator only compares the begin positions of the journal entries to determine their order.
         *
         * @retval std::weak_ordering The result of the comparison:
         * @rst
         * |   ``std::weak_ordering::less`` -- if the begin position of ``lhs`` is less than the begin position of ``rhs``.
         * |   ``std::weak_ordering::greater`` -- if the begin position of ``lhs`` is greater than the begin position of ``rhs``.
         * |   ``std::weak_ordering::equivalent`` -- if the begin position of ``lhs`` is equal to the begin position of ``rhs``.
         * @endrst
         */
        friend std::weak_ordering operator<=>(journal_entry const &lhs, journal_entry const &rhs) noexcept
        {
            return lhs.begin_position() <=> rhs.begin_position();
        }
        ///@}

        /** @name Member variables */
        ///@{
    private:
        /** @brief The begin position of the referenced segment in the journaled sequence. */
        size_type _position{};
        /** @brief The referenced segment in the journaled sequence. */
        sequence_type _segment{};
        ///@}
    };

    /** @name Deduction guides */
    ///@{
    /**
     * @relates journal_entry
     *
     * @tparam sequence_size_t The type of the position in the sequence. Must model:
     * \rst
     * |   ``std::integral``
     * \endrst
     * @tparam sequence_t The type of the sequence segment. Must model:
     * @rst
     * |   ``std::ranges::viewable_range``
     * |   ``std::ranges::random_access_range``
     * |   ``std::ranges::sized_range``
     * @endrst
     *
     * Deduces the template parameters from the given constructor arguments.
     */
    template <std::integral sequence_size_t, std::ranges::viewable_range sequence_t>
        requires std::ranges::random_access_range<sequence_t> &&
                 std::ranges::sized_range<sequence_t>
    journal_entry(sequence_size_t, sequence_t &&) -> journal_entry<sequence_size_t, std::ranges::iterator_t<sequence_t>>;
    ///@}

} // namespace libjst
