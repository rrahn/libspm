// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the node interface to be used with the lazy tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/utility/detail/multi_invocable.hpp>
#include <seqan3/range/views/type_reduce.hpp>

#include <libjst/set/concept_set.hpp>
#include <libjst/concept.hpp>
#include <libjst/journal.hpp>
#include <libjst/tree/jst_node_branch.hpp>

namespace libjst
{
    template <typename branch_state_t, typename variant_iterator_t>
    class jst_node_base
    {
    private:

        branch_state_t _state{}; // the one we need to create.
        variant_iterator_t _next_variant{};
        variant_iterator_t _last_variant{};
        // some internal values
        int32_t _max_branch_size{};
        int32_t _context_size{};

    public:
        // what kind of node type?
        using branch_node_type = jst_node_branch<branch_state_t, variant_iterator_t>;
        using value_type = branch_state_t;

        jst_node_base() = default;

        // how can we simplify this?
        explicit jst_node_base(branch_state_t state,
                               variant_iterator_t next_variant,
                               variant_iterator_t last_variant,
                               int32_t const max_branch_size,
                               size_t const context_size) :
            _state{std::move(state)},
            _next_variant{std::move(next_variant)},
            _last_variant{std::move(last_variant)},
            _max_branch_size{max_branch_size},
            _context_size{context_size - 1}
        {
            assert(context_size > 0);

            if (_next_variant != _last_variant)
                _state.set_range(0, libjst::position(_next_variant));  // current range to show.
            else  // so there must be an indicator when the state has no end anymore.
                _state.set_range(0, _max_branch_size);

            //     _next = libjst::position(*_next_variant);
            //     _last = _next + std::ranges::size(libjst::insertion(*_next_variant)) + _context_size;
            // } else {
            //     _next = std::ranges::size(_state.sequence());
            //     _last = _next;
            // }
        }

        value_type const & operator*() const noexcept
        {
            return _state;
        }

        value_type const * operator->() const noexcept
        {
            return std::addressof(_state);
        }

        bool has_next() const noexcept
        {
            return _next_variant != _last_variant;
        }

        // the increment function
        branch_node_type next()
        {
            assert(!has_next());

            // we are initialised to the beginning and have the first node
            auto branch_position = libjst::position(*_next_variant);
            auto branch_end_position = branch_position + std::ranges::ssize(libjst::insertion(*_next_variant)) +
                                        _context_size;
            branch_node_type variant_root{_state, _next_variant, _last_variant, branch_end_position};

            auto next_position = (++_next_variant != _last_variant) ? libjst::position(_next_variant) : _max_branch_size;
            assert(next_position >= branch_position);
            _state.set_range(branch_position, next_position);

            // move to the next
            // ++_next_variant; // move to the next variant
            // if (++_next_variant != _last_variant)
            //     _state.update(*_next_variant);  // or update, sufficient with one parameter?

            // if (_next_variant != _last_variant) {
            //     _next = libjst::position(*_next_variant);
            //     _last = _next + std::ranges::size(libjst::insertion(*_next_variant)) + _context_size;
            // } else {
            //     _next = std::ranges::size(_state.sequence());
            //     _last = std::ranges::size(_state.sequence());
            // }
            return variant_root;
        }

        constexpr operator bool() const noexcept
        {
            return hstate();
        }
    };
} // namespace libjst
