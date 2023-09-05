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

// using test_t = libjst::delta_sequence_variant<std::vector<seqan::alphabet_adaptor<seqan3::dna5> > >;
// using test_t2 = std::ranges::reverse_view<libjst::delta_sequence_variant<std::vector<seqan::alphabet_adaptor<seqan3::dna5> > > >;
// static_assert(std::ranges::viewable_range<test_t const &>);
// static_assert(std::ranges::viewable_range<test_t2 const &>);

namespace jstmap
{
    namespace detail {
        template <typename base_node_t>
        class unwind_node : public base_node_t {
        public:

            unwind_node(base_node_t && base) : base_node_t{std::move(base)}
            {}

            using base_node_t::reset_low;
        };
    } // namespace detail

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
        constexpr void operator()(cargo_t && base_cargo,
                                  finder_t && base_finder,
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

            // Seek to node whose label indicated a candidate region.
            auto seek_position = to_reverse_position(base_cargo.position());
            node_t initial_node = extend_tree.seek(seek_position);
            initial_node.toggle_alternate_path(); // Mark as alternate path to forbid extension on reference path

            // Set the correct right extension point, only because we know we left extended by this size
            // TODO: Do this on a level above
            // Alternative: Allow to get base journal sequence from the label and iterator so that we can get a position
            // on the root path which must be the same!
            // not the same sequences!
            auto reference_size = std::ranges::ssize(_reverse_tree.data().source());
            auto haystack_end = std::ranges::next(hostIterator(hostIterator(hostIterator(base_finder))), beginPosition(base_finder));
            auto prefix_begin_offset = std::ranges::distance(base_cargo.path_sequence().begin(), haystack_end);
            auto reverse_begin_offset = reference_size - 1 - prefix_begin_offset;

            auto initial_cargo = *initial_node;
            auto suffix_begin = std::ranges::next(initial_cargo.path_sequence().begin(), reverse_begin_offset);
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

                auto cargo = *node;
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

        constexpr libjst::seek_position to_forward_position(libjst::seek_position rev_seek_position) const noexcept {
            auto const & breakends = _reverse_tree.data().variants();
            return rev_seek_position.visit([&] (auto path_descriptor) {
                auto it = std::ranges::next(std::ranges::begin(breakends), rev_seek_position.get_variant_index());
                return unwind(std::move(path_descriptor), std::move(it));
            });
        }

    private:

        template <typename reverse_breakend_t>
        constexpr std::ptrdiff_t get_forward_index(reverse_breakend_t it) const noexcept {
            auto const & breakends = _reverse_tree.data().variants();
            return std::ranges::ssize(_reverse_tree.data().variants()) - std::ranges::distance(breakends.begin(), it);
        }

        constexpr libjst::seek_position to_reverse_position(libjst::seek_position fwd_seek_position) const noexcept {
            auto const & variants = _reverse_tree.data().variants();
            size_t breakend_count = variants.size();
            libjst::seek_position rev_seek_position{};
            fwd_seek_position.visit(seqan3::detail::multi_invocable{
                [&] (libjst::breakpoint_end) {
                    size_t right_breakend_idx = fwd_seek_position.get_variant_index();
                    libjst::breakpoint_end right_site = (*std::ranges::next(variants.begin(), right_breakend_idx)).get_breakpoint_end();
                    rev_seek_position.reset(breakend_count - right_breakend_idx, right_site);
                },
                [&] (libjst::alternate_path_descriptor const &) {
                    size_t right_breakend_idx = fwd_seek_position.get_variant_index();
                    rev_seek_position.initiate_alternate_node(breakend_count - right_breakend_idx - 1);
                }
            });
            return rev_seek_position;
        }

        template <typename breakend_iterator_t>
        constexpr libjst::seek_position unwind(libjst::breakpoint_end site, breakend_iterator_t it) const noexcept {
            libjst::seek_position fwd_position{};
            fwd_position.reset(get_forward_index(std::move(it)), site);
            return fwd_position;
        }

        template <typename breakend_iterator_t>
        constexpr libjst::seek_position unwind(libjst::alternate_path_descriptor const & descriptor,
                                               breakend_iterator_t it) const {
            auto unwind_tree = _reverse_tree | libjst::labelled() | libjst::merge();

            detail::unwind_node tmp{unwind_tree.root()};
            std::vector<std::ptrdiff_t> rev_path{};

            // Go to low breakpoint boundary
            rev_path.push_back(get_forward_index(it)); // record initial position
            --it;

            libjst::breakpoint_end low_end = (*it).get_breakpoint_end();
            tmp.reset_low(libjst::breakend_site<breakend_iterator_t>{std::move(it), low_end});
            for (auto it = std::ranges::begin(descriptor); it != std::ranges::end(descriptor); ++it) {
                if (*it) {
                    tmp = *tmp.next_alt();
                    rev_path.push_back(get_forward_index(tmp.low_boundary().get_breakend()));
                } else {
                    tmp = *tmp.next_ref();
                }
            }
            // Step 2: transform indices to fwd_position.
            auto fwd_path = rev_path | std::views::reverse;
            auto fwd_path_it = fwd_path.begin();
            auto fwd_path_end = fwd_path.end();

            std::ptrdiff_t last_index{*fwd_path_it};
            libjst::seek_position fwd_position{};
            fwd_position.initiate_alternate_node(last_index);
            ++fwd_path_it;

            // This is not true, because the offset is based on the nodes and not on the variants.
            for (; fwd_path_it != fwd_path_end; ++fwd_path_it) {
                std::ptrdiff_t jump_size = *fwd_path_it - last_index;
                // we must skip right breakends of the deletions!
                for (std::ptrdiff_t i = 0; i < jump_size; ++i) {
                    fwd_position.next_alternate_node(false);
                }
                fwd_position.next_alternate_node(true);
                last_index = *fwd_path_it;
            }

            return fwd_position;
        }
    };
}  // namespace jstmap
