// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides right extender.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <stack>

#include <libjst/matcher/myers_prefix_matcher_restorable.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

namespace jstmap
{
    template <typename base_tree_t, typename needle_t>
    class seed_suffix_extender
    {
        base_tree_t const & _base_tree;
        needle_t _needle;
        uint32_t _error_count{};
    public:
        seed_suffix_extender(base_tree_t const & base_tree, needle_t needle, uint32_t error_count) noexcept :
            _base_tree{base_tree},
            _needle{(needle_t &&) needle},
            _error_count{error_count}
        {}

        template <typename base_cargo_t, typename finder_t, typename callback_t>
        constexpr void operator()(base_cargo_t && base_cargo,
                                  finder_t && base_finder,
                                  callback_t && callback) const
        {
            libjst::restorable_myers_prefix_matcher extender{_needle, _error_count};

            auto extend_tree = _base_tree | libjst::labelled()
                                          | libjst::coloured()
                                          | libjst::trim(libjst::window_size(extender) - 1)
                                          | libjst::prune()
                                          | libjst::merge()
                                          | libjst::seek();

            using tree_t = decltype(extend_tree);
            using node_t = libjst::tree_node_t<tree_t>;
            using extender_t = decltype(extender);
            using cargo_t = libjst::tree_label_t<tree_t>;
            using state_t = typename extender_t::state_type;

            // Seek to node whose label indicated a candidate region.
            node_t initial_node = extend_tree.seek(base_cargo.position());
            initial_node.toggle_alternate_path(); // Mark as alternate path to forbid extension on reference path

            // Set the correct right extension point, only because we know we left extended by this size
            // TODO: Do this on a level above
            // Alternative: Allow to get base journal sequence from the label and iterator so that we can get a position
            // on the root path which must be the same!
            // not the same sequences!
            auto haystack_end = std::ranges::next(hostIterator(hostIterator(hostIterator(base_finder))), endPosition(base_finder));

            auto suffix_begin_offset = std::ranges::distance(base_cargo.path_sequence().begin(), haystack_end);
            auto initial_cargo = *initial_node;
            auto suffix_begin = std::ranges::next(initial_cargo.path_sequence().begin(), suffix_begin_offset);
            auto label_suffix = std::ranges::subrange{suffix_begin, initial_cargo.sequence().end()};

            // we call the extender and he updates its internal state!
            extender(label_suffix, [&] (auto && extender_finder) {
                callback(*initial_node, extender_finder, -getScore(extender.capture()));
            });

            // We do not need to extend the tree!
            if (std::ranges::size(label_suffix) >= libjst::window_size(extender) - 1)
                return;

            // Follow the path:
            std::stack<node_t> node_stack{};
            std::stack<state_t> state_stack{};

            // we may need some other comparison, because we don't want to move further.
            if (auto c_ref = initial_node.next_ref(); c_ref.has_value()) {
                node_stack.push(std::move(*c_ref));
                state_stack.push(extender.capture());
            }
            if (auto c_alt = initial_node.next_alt(); c_alt.has_value()) {
                node_stack.push(std::move(*c_alt));
                state_stack.push(extender.capture());
            }

            // Now we have to move here!
            while (!node_stack.empty()) {
                node_t node = std::move(node_stack.top());
                extender.restore(std::move(state_stack.top()));
                node_stack.pop();
                state_stack.pop();

                cargo_t cargo = *node;
                extender(cargo.sequence(), [&] (auto && extender_finder) {
                    callback(cargo, extender_finder, -getScore(extender.capture()));
                });

                if (auto c_ref = node.next_ref(); c_ref.has_value()) {
                    node_stack.push(std::move(*c_ref));
                    state_stack.push(extender.capture());
                }
                if (auto c_alt = node.next_alt(); c_alt.has_value()) {
                    node_stack.push(std::move(*c_alt));
                    state_stack.push(extender.capture());
                }
            }
        }
    };
}  // namespace jstmap
