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

#include <ranges>
#include <stack>

#include <libjst/matcher/myers_prefix_matcher_restorable.hpp>
#include <libjst/rcms/rcs_store_reversed.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

#include <jstmap/search/seed_prefix_node_cargo.hpp>
#include <jstmap/search/seed_prefix_seek_position.hpp>
#include <jstmap/search/seed_prefix_finder.hpp>

namespace jstmap
{
    template <typename base_tree_t, typename needle_t>
    class seed_prefix_extender
    {
        using reverse_rcs_t = decltype(libjst::rcs_store_reversed{std::declval<base_tree_t &&>().data().variants()});
        using tree_t = decltype(std::declval<reverse_rcs_t>() | libjst::make_volatile());
        using reverse_needle_t = decltype(std::declval<needle_t&&>() | std::views::reverse);

        tree_t _reverse_tree;
        reverse_needle_t _reverse_needle;
        uint32_t _error_count{};
    public:
        seed_prefix_extender(base_tree_t const & base_tree, needle_t needle, uint32_t error_count) noexcept :
            _reverse_tree{libjst::rcs_store_reversed{base_tree.data().variants()} | libjst::make_volatile()},
            _reverse_needle{(needle_t &&) needle | std::views::reverse},
            _error_count{error_count}
        {}

        template <typename cargo_t, typename finder_t, typename callback_t>
        constexpr void operator()(cargo_t && seed_cargo,
                                  finder_t && seed_finder,
                                  callback_t && callback) const
        {
            libjst::restorable_myers_prefix_matcher extender{_reverse_needle, _error_count};

            auto extend_tree = _reverse_tree | libjst::labelled()
                                             | libjst::coloured()
                                             | libjst::trim(libjst::window_size(extender) - 1)
                                             | libjst::prune()
                                             | libjst::merge()
                                             | libjst::seek();

            using tree_t = decltype(extend_tree);
            using node_t = libjst::tree_node_t<tree_t>;
            using extender_t = decltype(extender);
            using state_t = typename extender_t::state_type;

            auto const & variants = extend_tree.data().variants();
            seed_prefix_seek_position prefix_start_position{seed_cargo.position(), std::ranges::size(variants)};
            node_t initial_node = extend_tree.seek(std::move(prefix_start_position));
            initial_node.toggle_alternate_path(); // Mark as alternate path to forbid extension on reference path

            // Set the correct right extension point, only because we know we left extended by this size
            // TODO: Do this on a level above
            // Alternative: Allow to get base journal sequence from the label and iterator so that we can get a position
            // on the root path which must be the same!
            // not the same sequences!
            auto reference_size = std::ranges::ssize(_reverse_tree.data().source());
            auto haystack_end = std::ranges::next(hostIterator(hostIterator(hostIterator(seed_finder))), beginPosition(seed_finder));
            auto prefix_begin_offset = std::ranges::distance(seed_cargo.path_sequence().begin(), haystack_end);
            auto reverse_begin_offset = reference_size - 1 - prefix_begin_offset;

            seed_prefix_node_cargo initial_cargo{*initial_node, _reverse_tree};
            auto suffix_begin = std::ranges::next(initial_cargo.path_sequence().begin(), reverse_begin_offset);
            auto label_suffix = std::ranges::subrange{suffix_begin, initial_cargo.sequence().end()};

            // we call the extender and he updates its internal state!
            extender(label_suffix, [&] (auto const & extender_finder) {
                callback(initial_cargo, seed_prefix_finder{extender_finder, reference_size}, -getScore(extender.capture()));
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

                seed_prefix_node_cargo cargo{*node, _reverse_tree};
                extender(cargo.sequence(), [&] (auto const & extender_finder) {
                    callback(cargo, seed_prefix_finder{extender_finder, reference_size}, -getScore(extender.capture()));
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
