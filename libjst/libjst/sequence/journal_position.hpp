// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of journal position.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iterator>

namespace libjst
{

    template <typename journal_iterator_t>
    struct journal_position
    {
        // ----------------------------------------------------------------------------
        // Member Types
        // ----------------------------------------------------------------------------
    private:
        using journal_entry_type = std::iter_value_t<journal_iterator_t>;
        using sequence_type = typename journal_entry_type::sequence_type;

    public :

        using journal_iterator = journal_iterator_t;
        using sequence_iterator = std::ranges::iterator_t<sequence_type>;
        using size_type = typename journal_entry_type::size_type;

        // ----------------------------------------------------------------------------
        // Member Variables
        // ----------------------------------------------------------------------------

    public :

        journal_iterator journal_it{};
        sequence_iterator sequence_it{};

        // ----------------------------------------------------------------------------
        // Constructors, assignment and destructor

    public:
        constexpr journal_position() = default;
        constexpr journal_position(journal_iterator journal_it) :
            journal_position{journal_it, journal_it->segment().begin()}
        {
        }

        constexpr journal_position(journal_iterator journal_it, sequence_iterator sequence_it) :
            journal_it{std::move(journal_it)},
            sequence_it{std::move(sequence_it)}
        {
        }

        // ----------------------------------------------------------------------------
        // Utility functions
    private:
        friend bool operator==(journal_position const &lhs, journal_position const &rhs) noexcept = default;
        friend auto operator<=>(journal_position const &lhs, journal_position const &rhs) noexcept = default;

        friend size_type to_sequence_position(journal_position const &pos) noexcept
        {
            auto [journal_it, sequence_it] = pos;
            return journal_it->begin_position() + std::ranges::distance(journal_it->segment().begin(), sequence_it);
        }

        friend auto split_at(journal_position const &pos) noexcept
        {
            return pos.journal_it->split_at(pos.sequence_it);
        }
    };
} // namespace libjst
