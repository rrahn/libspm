// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides breakpoint node.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/sequence_tree/breakend_site.hpp>
#include <libjst/sequence_tree/node_descriptor.hpp>

namespace libjst
{
    template <typename breakend_iterator>
    class breakpoint_node : public node_descriptor
    {
    private:
        using position_t = breakend_site<breakend_iterator>;

        position_t _low{};
        position_t _high{};

    public:
        using position_type = position_t;
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr breakpoint_node() = default; //!< Default.
        explicit constexpr breakpoint_node(position_type low, position_type high) :
            _low{std::move(low)},
            _high{std::move(high)}
        {
        }
        //!\}

        constexpr position_t const & low_boundary() const noexcept {
            return _low;
        }

        constexpr position_t const & high_boundary() const noexcept {
            return _high;
        }

        constexpr std::optional<breakpoint_node> next_alt() const noexcept {
            if (this->from_reference() && high_boundary().is_low_end()) {
                position_t child_low = high_boundary();
                position_t child_high = next_high_boundary_alt(high_boundary());
                breakpoint_node child{std::move(child_low), std::move(child_high)};
                child.activate_state(node_state::variant);
                return child;
            } else {
                return std::nullopt;
            }
        }

        constexpr breakpoint_node next_ref() const noexcept {
            // either from alternate path
            position_t child_low = high_boundary();
            position_t child_high = next_high_boundary_ref(high_boundary());
            breakpoint_node child{std::move(child_low), std::move(child_high)};
            if (this->on_alternate_path())
                child.toggle_alternate_path();
            return child;

            // if (node_descriptor::from_reference()) {
            //     if (is_leaf_of_alternate_subtree()) {
            //         return std::nullopt;
            //     } else {
            //         return visit_next_ref_impl(std::move(child));
            //     }
            // } else {
            //     assert(get_left() == get_right());
            //     child.set_left(get_right());
            //     child.set_right(find_next_valid_right_variant());
            //     child.set_next(next_variant_after(child.get_right()));
            //     child.initialise_reference_state_from();
            // }
            // return child;
        }

    private:

        position_t next_high_boundary_alt(position_t const & boundary) const noexcept {
            if (boundary.is_low_end()) {
                auto breakend = (*boundary).jump_to_mate().value_or(boundary.get_breakend());
                return position_t{std::move(breakend), breakpoint_end::high};
            }
            return next_high_boundary_ref(boundary);
        }

        position_t next_high_boundary_ref(position_t const & boundary) const noexcept {
            auto next_breakend = std::ranges::next(boundary.get_breakend());
            breakpoint_end next_site = (*next_breakend).get_breakpoint_end();
            position_t next_high_boundary{std::move(next_breakend), std::move(next_site)};
            // Move over all
            while (libjst::position(high_boundary()) > libjst::position(next_high_boundary)) {
                next_breakend = std::ranges::next(next_high_boundary.get_breakend());
                next_site = (*next_breakend).get_breakpoint_end();
                next_high_boundary = position_t{std::move(next_breakend), std::move(next_site)};
            }
            return next_high_boundary;
        }

        friend constexpr bool operator==(breakpoint_node const &, breakpoint_node const &) noexcept = default;
    };
}  // namespace libjst
