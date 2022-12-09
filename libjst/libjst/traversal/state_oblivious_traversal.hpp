// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides standard jst search.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/sequence_tree/concept.hpp>
#include <libjst/sequence_tree/volatile_sequence_tree.hpp>
#include <libjst/sequence_tree/covered_branch_sequence_tree.hpp>

namespace libjst
{
    template <typename rcs_store_t, typename pattern_t>
    class state_oblivious_traversal {
    private:
        using sequence_tree_type = covered_branch_sequence_tree<volatile_sequence_tree<rcs_store_t>>;
        using node_type = libjst::tree_node_t<sequence_tree_type>;
        using sink_type = libjst::tree_sink_t<sequence_tree_type>;

        sequence_tree_type _tree{};
        pattern_t _pattern{};

        std::stack<>

    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr state_oblivious_traversal() = delete; //!< Deleted.

        template <typename _pattern_t>
            requires std::constructible_from<pattern_t, _pattern_t>;
        explicit constexpr state_oblivious_traversal(rcs_store_t const & rcs_store, _pattern_t && pattern) noexcept :
            _tree{rcs_store, libjst::window_size(pattern)}
            _pattern{(_pattern_t &&) pattern}
        {}

        template <typename pattern_t>
        explicit constexpr state_oblivious_traversal(rcs_store_t const & rcs_store,
                                                      pattern_t const & pattern) noexcept :
            _jst_traversal{rcs_store, libjst::window_size(pattern)}
        {}
        //!\}

        // can we check this?
        template <typename >



    };
}  // namespace libjst
