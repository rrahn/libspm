// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::reference_position.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>

namespace libjst
{
/*!\brief A specific position type for the libjst::journaled_sequence_tree reference collecetion.
 *
 * \details
 *
 * This position can be used to index any position within the underlying reference collection.
 * The index identifies the reference within the set and the offset gives the position whithin this
 * reference sequence.
 */
struct reference_position
{
    size_t idx{}; //!< The index of the reference sequence.
    size_t offset{}; //!< The offset within the respective reference sequence.

    reference_position & operator++() noexcept
    {
        ++offset;
        return *this;
    }

    reference_position & operator+=(size_t const n) noexcept
    {
        offset += n;
        return *this;
    }

    reference_position operator+(size_t const n) const noexcept
    {
        reference_position tmp{*this};
        tmp.offset += n;
        return tmp;
    }

    friend reference_position operator+(size_t const n, reference_position const & rhs) noexcept
    {
        return rhs + n;
    }

    reference_position & operator--() noexcept
    {
        --offset;
        return *this;
    }

    reference_position & operator-=(size_t const n) noexcept
    {
        offset -= n;
        return *this;
    }

    reference_position operator-(size_t const n) const noexcept
    {
        reference_position tmp{*this};
        tmp.offset -= n;
        return tmp;
    }

    //!\brief Default three-way comparison.
    constexpr auto operator<=>(reference_position const &) const noexcept = default;
};

/*!\brief A specific position type for the libjst::journaled_sequence_tree reference collection.
 * \relates libjst::reference_position
 *
 * \tparam stream_t The type of the output stream.
 * \tparam reference_position_t The type of the reference postion; must be the same as libjst::reference_position.
 *
 * \param[in] stream The stream to write to.
 * \param[in] reference_pos The object to stream.
 */
template <typename char_t, typename char_traits_t, typename reference_position_t>
//!\cond
    requires std::same_as<std::remove_cvref_t<reference_position_t>, reference_position>
//!\endcond
inline std::basic_ostream<char_t, char_traits_t> & operator<<(std::basic_ostream<char_t, char_traits_t> & stream, reference_position_t && reference_pos)
{
    stream << "["
                << "idx: " << reference_pos.idx << ", "
                << "pos: " << reference_pos.offset
           << "]";
    return stream;
}
}  // namespace libjst
