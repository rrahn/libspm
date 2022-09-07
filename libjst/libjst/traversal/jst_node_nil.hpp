// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides nil node for the jst.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace libjst
{

    template <typename linked_node_t>
    struct jst_node_nil
    {
        constexpr bool operator==(linked_node_t const & node) const noexcept
        {
            // I need some kind of signal to now when the tree ends!
            return node.at_end();
        }
    };
}  // namespace libjst
