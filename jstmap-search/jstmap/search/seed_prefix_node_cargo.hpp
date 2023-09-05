// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seek position.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>
#include <libjst/sequence_tree/seek_position.hpp>

namespace jstmap
{
    namespace detail {
        template <typename base_node_t>
        class unwind_node : public base_node_t {
        public:

            unwind_node(base_node_t && base) : base_node_t{std::move(base)}
            {}

            using base_node_t::reset;
        };
    } // namespace detail

    template <typename cargo_t, typename reverse_tree_t>
    class seed_prefix_node_cargo : public cargo_t {
    private:
        using base_t = cargo_t;

        reverse_tree_t const & _reverse_tree;

    public:

        explicit constexpr seed_prefix_node_cargo() = delete;
        explicit constexpr seed_prefix_node_cargo(cargo_t cargo, reverse_tree_t const & reverse_tree) noexcept :
            base_t{std::move(cargo)},
            _reverse_tree{reverse_tree}
        {}

        constexpr libjst::seek_position position() const noexcept {
            return to_forward_position(base_t::position());
        }

    private:

        constexpr libjst::seek_position to_forward_position(libjst::seek_position reverse_position) const noexcept {
            auto const & breakends = _reverse_tree.data().variants();
            return reverse_position.visit([&] (auto path_descriptor) {
                auto it = std::ranges::next(std::ranges::begin(breakends), reverse_position.get_variant_index());
                return unwind(std::move(path_descriptor), std::move(it));
            });
        }

        template <typename reverse_breakend_t>
        constexpr std::ptrdiff_t get_forward_index(reverse_breakend_t it) const noexcept {
            auto const & breakends = _reverse_tree.data().variants();
            return std::ranges::ssize(_reverse_tree.data().variants()) - std::ranges::distance(breakends.begin(), it) - 1;
        }

        template <typename breakend_iterator_t>
        constexpr libjst::seek_position unwind(libjst::breakpoint_end site, breakend_iterator_t it) const noexcept {
            libjst::seek_position fwd_position{};
            fwd_position.reset(get_forward_index(std::move(it)) - 1, site);
            return fwd_position;
        }

        template <typename breakend_iterator_t>
        constexpr libjst::seek_position unwind(libjst::alternate_path_descriptor const & descriptor,
                                               breakend_iterator_t it) const {
            auto unwind_tree = _reverse_tree | libjst::labelled() | libjst::merge() | libjst::seek();

            detail::unwind_node tmp{unwind_tree.root()};
            std::vector<std::ptrdiff_t> rev_path{};

            // Go to low breakpoint boundary
            // rev_path.push_back(get_forward_index(it)); // record initial position
            --it;
            libjst::breakpoint_end low_end = (*it).get_breakpoint_end();
            libjst::seek_position initial_position{};
            initial_position.reset(get_forward_index(it), low_end);
            tmp.reset(libjst::breakend_site<breakend_iterator_t>{std::move(it), low_end}, std::move(initial_position));
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
                assert(*fwd_path_it > last_index);
                std::ptrdiff_t skipped_breakends = *fwd_path_it - last_index - 1;
                // we must skip right breakends of the deletions!
                for (std::ptrdiff_t i = 0; i < skipped_breakends; ++i) {
                    fwd_position.next_alternate_node(false);
                }
                fwd_position.next_alternate_node(true);
                last_index = *fwd_path_it;
            }

            return fwd_position;
        }
    };
}  // namespace jstmap
