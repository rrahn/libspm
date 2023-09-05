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
#include <libjst/sequence_tree/partial_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

#include <jstmap/global/match_position.hpp>
#include <jstmap/search/seed_node_wrapper.hpp>
#include <jstmap/search/seed_extension_tree.hpp>

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

        template <typename seed_cargo_t, typename finder_t, typename callback_t>
        constexpr void operator()(seed_cargo_t && seed_cargo,
                                  finder_t && seed_finder,
                                  callback_t && callback) const
        {
            if (std::ranges::empty(_needle)) {
                callback(match_position{.tree_position = seed_cargo.position(),
                                        .label_offset = endPosition(seed_finder)},
                                        _error_count);
                return;
            }
            libjst::restorable_myers_prefix_matcher extender{_needle, _error_count};

            // A) Is inside an alternate path -> we can set the partial end to the variant index
            //                                -> that is what we extended it with
            // But when we are here, then we are already inside a path and we can not go back into the tree
                // base node is already a natural end -> so we cannot go back.
                // - So there must be something else?

            std::ptrdiff_t distance_to_end = std::ranges::ssize(seed_cargo.sequence()) - endPosition(seed_finder);
            match_position start{.tree_position = seed_cargo.position(),
                                 .label_offset = std::ranges::ssize(seed_cargo.path_sequence()) - distance_to_end};

            auto tree_adapter = libjst::labelled()
                              | libjst::coloured()
                              | libjst::prune()
                              | libjst::merge()
                              | jstmap::extend_from(start, libjst::window_size(extender) - 1)
                              | libjst::seek();

            auto extend_tree = _base_tree | tree_adapter;
            // auto [on_alternate_path, extend_tree] = seed_cargo.position().visit(seqan3::detail::multi_invocable{
            //     [&] (libjst::breakpoint_end site) {
            //         auto const & rcms = _base_tree.data();
            //         return std::pair{false, libjst::partial_tree(rcms, endPosition(seed_finder), 0) | tree_adapter()};
            //     },
            //     [&] (...) {
            //         return std::pair{true, _base_tree | tree_adapter()};
            //     }
            // });

            // using tree_t = decltype(extend_tree);
            // using node_t = libjst::tree_node_t<tree_t>;
            // using extender_t = decltype(extender);
            // using cargo_t = libjst::tree_label_t<tree_t>;
            // using state_t = typename extender_t::state_type;

            // Seek to node whose label indicated a candidate region.
            // we may only seek in case of a alternate node!
            // node_t initial_node{};
            // if (on_alternate_path)
            //     initial_node = extend_tree.seek(seed_cargo.position());
            // else
            //     initial_node = extend_tree.root();

            // Case a: inital_node = partial root
                // - begin position already set
                // - extension is set as well
                // - run into last alternate but not reference
                // - do normal traversal
            // Case b: initial_node = rewind path but position not set
                // - a new sequence label, with new start offset?
                // - is in the last node of alternate path with maybe empty suffix.
                // -
            // if (on_alternate_path) {
            //     // so we get the haystack position
            //     // what we want to know, how far is this end position away from the root.
            //     // because when going to the same node we must be able to set to the same position
            //     auto suffix_start_it = hostIterator(hostIterator(end(seed_finder))); // iterator to local haystack segment of current node
            //     std::ptrdiff_t suffix_start_position = std::ranges::distance(seed_cargo.path_sequence().begin(), suffix_start_it);

            //     auto initial_cargo = *initial_node;
            //     auto suffix_begin = std::ranges::next(initial_cargo.path_sequence().begin(), suffix_start_position);
            //     auto label_suffix = std::ranges::subrange{suffix_begin, initial_cargo.path_sequence().end()};

            //     // we call the extender and he updates its internal state!
            //     extender(label_suffix, [&] (auto && suffix_finder) {
            //         callback(match_position{.tree_position = initial_cargo.position(),
            //                                 .label_offset = endPosition(suffix_finder)},
            //                                 -getScore(extender.capture()));
            //     });
            // }

            // We do not need to extend the tree!
            // if (std::ranges::size(label_suffix) >= libjst::window_size(extender) - 1)
            //     return;

            libjst::tree_traverser_base suffix_traverser{extend_tree};
            extension_state_manager manager{extender};
            suffix_traverser.subscribe(manager);
            for (auto cargo : suffix_traverser) {
                // auto cargo = *node;
                extender(cargo.sequence(), [&] (auto && suffix_finder) {
                    callback(match_position{.tree_position = cargo.position(),
                                            .label_offset = endPosition(suffix_finder)},
                                            -getScore(extender.capture()));
                });
            }

            // // Follow the path:
            // std::stack<node_t> node_stack{};
            // std::stack<state_t> state_stack{};

            // // we may need some other comparison, because we don't want to move further.
            // if (auto c_ref = initial_node.next_ref(); c_ref.has_value()) {
            //     node_stack.push(std::move(*c_ref));
            //     state_stack.push(extender.capture());
            // }
            // if (auto c_alt = initial_node.next_alt(); c_alt.has_value()) {
            //     node_stack.push(std::move(*c_alt));
            //     state_stack.push(extender.capture());
            // }

            // // Now we have to move here!
            // while (!node_stack.empty()) {
            //     node_t node = std::move(node_stack.top());
            //     extender.restore(std::move(state_stack.top()));
            //     node_stack.pop();
            //     state_stack.pop();

            //     cargo_t cargo = *node;
            //     extender(cargo.sequence(), [&] (auto && suffix_finder) {
            //         callback(match_position{.tree_position = cargo.position(),
            //                                 .label_offset = endPosition(suffix_finder)},
            //                                 -getScore(extender.capture()));
            //     });

            //     if (auto c_ref = node.next_ref(); c_ref.has_value()) {
            //         node_stack.push(std::move(*c_ref));
            //         state_stack.push(extender.capture());
            //     }
            //     if (auto c_alt = node.next_alt(); c_alt.has_value()) {
            //         node_stack.push(std::move(*c_alt));
            //         state_stack.push(extender.capture());
            //     }
            // }
        }
    };
}  // namespace jstmap
