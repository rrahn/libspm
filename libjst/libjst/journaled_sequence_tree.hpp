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

#include <seqan3/std/algorithm>
#include <seqan3/std/filesystem>
#include <iostream>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

#include <cereal/types/vector.hpp>

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/range/concept.hpp>

namespace libjst::no_adl
{
//!\brief Specific class implementation in no_adl namespace to avoid ADL of template arguments.
template <seqan3::sequence sequence_t>
class journaled_sequence_tree
{
public:
    class type;
};

//!\brief Implements the actual journaled sequence tree type.
template <seqan3::sequence sequence_t>
class journaled_sequence_tree<sequence_t>::type
{
private:
    sequence_t _reference; //!< The internal reference used for referential compression.
    std::vector<sequence_t> _sequences; //!< The stored sequences.

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr type() = default; //!< Default.
    constexpr type(type const &) = default; //!< Default.
    constexpr type(type &&) = default; //!< Default.
    constexpr type & operator=(type const &) = default; //!< Default.
    constexpr type & operator=(type &&) = default; //!< Default.
    ~type() = default; //!< Default.

    /*!\brief Constructs the journaled sequence tree with a given reference sequence.
     *
     * \param[in] reference The reference sequence to set.
     *
     * \details
     *
     * The journaled sequence tree takes ownership over the passed reference sequence. Accordingly, only temporaries
     * or moved sequences can be used. To access the reference later one can use the
     * libjst::journaled_sequence_tree::reference member function to access the stored reference.
     * Thus, one can use this reference sequence as long as the journaled sequence tree is valid.
     */
    type(sequence_t && reference) : _reference{std::move(reference)}
    {}
    //!\}

    //!\brief Returns the stored reference sequence.
    sequence_t const & reference() const
    {
        return _reference;
    }

    //!\brief Returns the number of stored sequences.
    size_t size() const
    {
        return _sequences.size();
    }

    /*!\brief Adds a new sequence to the journaled sequence tree based on the given pairwise alignment.
     *
     * \tparam alignment_t The type of the alignment to add; must model pairwise_alignment.
     *
     * \param[in] alignment The alignment to add to this object.
     *
     * \details
     *
     * This method adds a sequence to the current journaled sequence tree using referential compression. The given
     * alignment is used to determine the transcript to derive the added sequence from the stored reference sequence.
     * Only pairwise alignments are supported.
     * The first sequence of the alignment must be identical to the sequence returned by
     * libjst::journaled_sequence_tree::reference after all gap characters have been removed.
     * The second sequence is added encoded by the given alignment.
     */
    template <typename alignment_t>
    void add(alignment_t && alignment)
    {
        auto & [ref, target] = alignment;

        constexpr auto without_gap  = [] (char const c) -> bool { return c != '-'; };

        auto pure_target_view = target | std::views::filter(without_gap);

        if (!std::ranges::equal(ref | std::views::filter(without_gap), _reference))
            throw std::invalid_argument{"The first aligned sequence must be equal to the reference sequence of this "
                                        "journaled sequence tree without the gaps."};

        sequence_t pure_target_sequence;
        std::ranges::copy(target | std::views::filter(without_gap), std::cpp20::back_inserter(pure_target_sequence));
        _sequences.push_back(std::move(pure_target_sequence));
    }

    /*!\brief Saves the journaled sequence tree to the given output archive.
     *
     * \tparam output_archive_t The type of the output_archive; must model seqan3::cereal_output_archive.
     *
     * \param[in] archive The archive to serialise this object to.
     */
    template <seqan3::cereal_output_archive output_archive_t>
    void save(output_archive_t & archive) const
    {
        archive(_reference, _sequences);
    }
};
} // namespace libjst::no_adl

namespace libjst
{

/*!\brief A referentially compressed sequence tree over collection of sequences.
 *
 * \tparam sequence_t The type of the sequences to store; must model seqan3::sequence.
 *
 * \details
 *
 * This class stores a collection of sequences in a referentially compressed way to tremendously reduce the
 * memory footprint for storing large collections of sequences with a high similarity. Sequences can be added by
 * adding an alignment between the stored reference sequence and the respective target sequence. This class further
 * supports a special cursor to enable an efficient, compression parallel traversal over the stored sequences.
 * This can be used in conjunction with any context based streaming algorithm to speed-up the search against large
 * collection of sequences.
 */
template <seqan3::sequence sequence_t>
using journaled_sequence_tree = typename no_adl::journaled_sequence_tree<sequence_t>::type;
}  // namespace libjst
