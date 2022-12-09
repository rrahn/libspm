// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides jst traversal with state obliviousness.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace libjst
{
    // TODO: require rcs_store interface.
    template <typename rcs_store_t>
    class jst_traverser_state_aggregation
    {
    private:

        // properties of the journaled sequence tree are coming from the rcs_store_t
        using jst_t = journaled_sequence_tree<rcs_store_t>;
    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr jst_traverser_state_aggregation() = default; //!< Default.

        explicit constexpr jst_traverser_state_aggregation(rcs_store_t const & rcs_store)
        {}
        //!\}





    };
}  // namespace libjst
