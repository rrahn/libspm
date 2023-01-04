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
#include <libjst/sequence_tree/left_extend_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/prune_unsupported.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

namespace libjst
{
    struct state_oblivious_traverser {
        template <typename tree_t, typename pattern_t, typename callback_t>
        constexpr void operator()(tree_t && tree, pattern_t && pattern, callback_t && callback) const {
            if (libjst::window_size(pattern) == 0)
                return;

            auto search_tree = tree | libjst::labelled<libjst::sequence_label_kind::root_path>()
                                               | libjst::coloured()
                                               | trim(libjst::window_size(pattern) - 1)
                                               | prune_unsupported()
                                               | left_extend(libjst::window_size(pattern) - 1)
                                               | merge(); // make big nodes

            tree_traverser_base oblivious_path{search_tree};
            for (auto it = oblivious_path.begin(); it != oblivious_path.end(); ++it) {
                auto && label = *it;
                pattern(label.sequence(), [&] (auto && label_it) {
                    callback(std::move(label_it), label); // either cargo offers access to node context or not!
                });
            }
        }
    };
}  // namespace libjst
