// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides referentially compressed sequence store searcher.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace libjst
{

    // TODO: Add constraints to rcs_store_t and search_strategy_t
    template <typename rcs_store_t, typename search_strategy_t>
    class rcs_store_searcher
    {
    private:

        rcs_store_t const & _rcs_store;
    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr rcs_store_searcher() = delete; //!< Deleted.
        explicit constexpr rcs_store_searcher(rcs_store_t const & rcs_store) noexcept : _rcs_store{rcs_store}
        {}
        //!\}

        template <typename pattern_t, typename callback_t>
        constexpr void operator()(pattern_t && pattern, callback_t && callback) const {
            // we need to choose the search strategy based on the pattern.

            // does the traverser have a state and return a location?
            // it would allow us to interrupt the search and not explore everything!
            // is like the finder concept, but smarter!
            // we may even use it to adapt the finder object from SeqAn?
            // ideally we like to get some type erased version here and can assign a concrete traverser instance, allocated in here.
            auto rcs_haystack = make_rcs_traverser<pattern_t>(_rcs_store); // only based on features of pattern?
            // let's say the traverser is a burrowed range?
            // who allocates the traversal?
            // We may customise the traverser by augmenting it with a special sentinel symbol.

            // making traverser a range interface over the haystack.
            // we can always set the pattern state as well, we need to get the default pattern state the first time
            // we call the operation.

            // the point is I can add an artificial end by selecting a position inside the tree that we can not exceed in the base sequence.
            // borrowed_subrange_t<traverser &> subrange_traverser{traverser};
            while (true) {
                if (rcs_haystack = libjst::search(rcs_haystack, pattern); !std::ranges::empty(rcs_haystack))
                    break;

                traverser.coordinate() // etc.
                report(std::move(it)); // maybe we need some extra infromation, that can only extracted from the traverser?                }
            }

            // We can orchestrate the partial jst search from outside!
            // that is we can make sure if we are passing some specific information about the tree.
            // jst is build over the full tree.

            // rules for searching which can be adapted through CPOs?
        }

    private:

        template <typename pattern_t>
            // requires ... pattern
        constexpr auto make_jst_traversal(pattern_t && pattern) const noexcept {
            return jst_traverser_state_obliviousness{_rcs_store, (pattern_t &&) pattern};
        }
    };
}  // namespace libjst
