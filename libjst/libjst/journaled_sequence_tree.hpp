// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::journaled_sequence_tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iostream>

#include <seqan3/std/filesystem>

namespace libjst
{
template <typename sequences_t>
class journaled_sequence_tree
{
private:

    sequences_t sequences{};
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr journaled_sequence_tree() = default; //!< Default.
    constexpr journaled_sequence_tree(journaled_sequence_tree const &) = default; //!< Default.
    constexpr journaled_sequence_tree(journaled_sequence_tree &&) = default; //!< Default.
    constexpr journaled_sequence_tree & operator=(journaled_sequence_tree const &) = default; //!< Default.
    constexpr journaled_sequence_tree & operator=(journaled_sequence_tree &&) = default; //!< Default.
    ~journaled_sequence_tree() = default; //!< Default.

    journaled_sequence_tree(sequences_t sequences) : sequences{std::move(sequences)}
    {}
    //!\}

    void save(std::filesystem::path const &)
    {
        std::cout << "Write me!\n";
    }
};
}  // namespace libjst
