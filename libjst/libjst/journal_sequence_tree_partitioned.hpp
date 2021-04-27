// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::journaled_sequence_tree partitioned.
 * \author Tom Lukas Lankenau <tom.lankenau AT fu-berlin.de>
 */

#pragma once

#include <math.h> // ceil
#include <vector> // vector

#include <libjst/journal_sequence_tree_context_enumerator.hpp>
#include <libjst/journaled_sequence_tree.hpp>

namespace libjst
{

template <typename jst_t>
class journal_sequence_tree_partitioned
{
private:
    /* data */
    using enumerator_type = typename jst_t::context_enumerator_type;

    std::vector<enumerator_type> _bins;
    const jst_t * _jst;
public:

    constexpr journal_sequence_tree_partitioned() = default; //!< Default.
    constexpr journal_sequence_tree_partitioned(journal_sequence_tree_partitioned const &) = default;
        //!< Default.
    constexpr journal_sequence_tree_partitioned(journal_sequence_tree_partitioned &&) = default;
        //!< Default.
    constexpr journal_sequence_tree_partitioned & operator=(journal_sequence_tree_partitioned const &)
        = default; //!< Default.
    constexpr journal_sequence_tree_partitioned & operator=(journal_sequence_tree_partitioned &&)
        = default; //!< Default.
    ~journal_sequence_tree_partitioned() = default; //!< Default.

    journal_sequence_tree_partitioned (const jst_t * jst, size_t bin_count):
    _bins(),
    _jst(jst)
    {
        assert(_jst != nullptr);

        size_t bin_size = ceil(_jst->reference().size() / bin_count);
        for (size_t i = 0; i < bin_count; i++)
        {
            _bins.push_back(enumerator_type(_jst, 1, i * bin_size, (i + 1) * bin_size));
        }

    }

    enumerator_type operator[](size_t i)
    {
        assert(i < _bins.size());
        return _bins[i];
    }

    size_t size()
    {
        return _bins.size();
    }
};

} // namespace libjst

// // member
// jst
// bin_count (?)
// bin_size  (?)
// jst_enumerator
//
// // member functions
// [i].begin() -> jst_enumerator[i]
// [i].end() -> jst_enumerator[i+1] bzw. jst.end();
// serialise() -> archive()
// }

// /*!\brief Saves this journaled sequence tree to the given output archive.
//  *
//  * \tparam output_archive_t The type of the output_archive; must model seqan3::cereal_output_archive.
//  *
//  * \param[in, out] archive The archive to serialise this object to.
//  */
// template <seqan3::cereal_output_archive output_archive_t>
// void save(output_archive_t & archive) const
// {
//     archive(_reference, _delta_events, _size);
// }
