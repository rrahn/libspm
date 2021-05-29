// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::journal_sequence_tree_coordinate.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace libjst
{

/*!\brief Represents a specific position within the journal sequence tree.
 *
 * \details
 *
 * Using this coordinate any position within the journal sequence tree can be efficienly recovered through a seek
 * operation. Jumps to the reference position and only the subtree needs to be parsed if applicable.
 */
struct journal_sequence_tree_coordinate
{
    uint32_t position{}; //!< The offset of the context head from the respective subtree root.
    uint32_t event_id{}; //!< The id of the next event.
    uint32_t subtree_steps{}; //!< Number of steps we have to advance into the subtree to get to the correct position.
    uint32_t context_size{}; //!< The context size used to step into the tree.
};

template <typename char_t, typename char_traits_t>
inline std::basic_ostream<char_t, char_traits_t> & operator<<(std::basic_ostream<char_t, char_traits_t> & stream,
                                                              journal_sequence_tree_coordinate const & coordinate)
{
    stream << "[" << coordinate.position <<
              ", " << coordinate.event_id <<
              ", " << coordinate.subtree_steps <<
              ", " << coordinate.context_size << "]";
    return stream;
}

}  // namespace libjst
