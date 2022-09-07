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

#include <concepts>

namespace libjst
{
    // assuming this one is burrowed?
    // then we need to get everything from the iterator.
    // but this means copying a stack or something that we don't want.
    // that is we can return a new range from the burrowed range here.
    // iterator could keep shared ptr of tree
    // traversal_context: includes the actual tree.
    // then we have encapsulated the tree operations from the algorithms
    // so one part could be search algorithm using online pattern
    // other could be extension algorithm
    // one could be different stuff as well.
    template <typename rcs_store_t>
    class jst_traverser_state_stashing // realising burrowed_range?
    {
    private:

        traversal_context_t _traversal_context{};
        state_stash<state_t> _stash{};

        class iterator;
        class sentinel;

    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr jst_traverser_state_stashing() = default; //!< Default.

        template <typename pattern_t>
        explicit constexpr jst_traverser_state_stashing(rcs_store_t const & rcs_store, pattern_t const & pattern)
        {
            // create jst
            // initial state from the pattern!
            _stash.cache(pattern.capture()); // initial snapshot captured in stash.
            // register _stash to jst notifier
        }
        //!\}

        constexpr iterator begin() noexcept {
            return iterator{}; // handle to traverser using move only operations! nothing to worry about?
        }

        constexpr sentinel end() noexcept {
            return sentinel{};
        }

    private:

        template <typename pattern_t>
            requires std::predicate<pattern_t, typename context_traits<traversal_context_t>::sequence_type>
        constexpr friend iterator tag_invoke(std::tag_t<search>, iterator it, sentinel const last, pattern_t && pattern) {
            for (; it != last; ++it) {
                pattern.restore(it->state()); // now we set the internal state to the pattern.
                if (std::invoke(pattern, it->sequence()))
                    break;
                it->state(pattern.capture());
            }
            it->state(pattern.capture());
            return it;
        }
    };
}  // namespace libjst
