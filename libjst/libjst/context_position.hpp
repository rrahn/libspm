// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::context_position.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>

namespace libjst
{
/*!\brief A specific position type for the libjst::journaled_sequence_tree_cursor.
 *
 * \details
 *
 * Provides a description of the context position of a particular sequence by specifying the sequence id and the
 * begin poisition of the respective sequence context within this sequence.
 */
struct context_position
{
    size_t sequence_id; //!< The id of the sequence.
    size_t sequence_position; //!< The context begin position within the sequence.

    //!\brief Default three-way comparison.
    constexpr auto operator<=>(context_position const &) const noexcept = default;
};

/*!\brief A specific position type for the libjst::journaled_sequence_tree_cursor.
 * \relates libjst::context_position
 *
 * \tparam stream_t The type of the output stream.
 * \tparam context_position_t The type of the context postion; must be the same as libjst::context_position.
 *
 * \param[in] stream The stream to write to.
 * \param[in] context_pos The object to serialise.
 */
template <typename stream_t, std::same_as<context_position> context_position_t>
stream_t & operator<<(stream_t & stream, context_position_t context_pos)
{
    stream << "["
                << "id: " << context_pos.sequence_id << ", "
                << "pos: " << context_pos.sequence_position
           << "]";
    return stream;
}
}  // namespace libjst
