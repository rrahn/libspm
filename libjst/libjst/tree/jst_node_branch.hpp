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

namespace libjst
{
    // Clearly needs some refactoring!
    template <typename branch_state_t, typename variant_iterator_t>
    class jst_node_branch
    {
    private:
        branch_state_t _state{};
        variant_iterator_t _next_variant{};
        variant_iterator_t _last_variant{};
        int32_t _max_branch_size{};
        int32_t _offset{};

    public:

        using value_type = branch_state_t;

        jst_node_branch() = default;
        jst_node_branch(jst_node_branch const &) = default;
        jst_node_branch(jst_node_branch &&) = default;
        jst_node_branch & operator=(jst_node_branch const &) = default;
        jst_node_branch & operator=(jst_node_branch &&) = default;

        // construction from base node
        explicit jst_node_branch(branch_state_t state, // parent_state
                                 variant_iterator_t const current_variant,
                                 variant_iterator_t last_variant,
                                 size_t const max_branch_size,
                                 int32_t const offset = 0) :
            _state{std::move(state)},
            _next_variant{current_variant},
            _last_variant{std::move(last_variant)},
            _max_branch_size{max_branch_size},
            _offset{offset}
        {
            assert(_next_variant != _last_variant);

            _state.set_branch(*_next_variant, _max_branch_size);
            auto const branch_position = libjst::position(*current_variant);
            auto const branch_end = branch_position + libjst::deletion(*current_variant);

            // First find first variant that is not an insertion including the current variant.
            _next_variant = std::ranges::find_if(++_next_variant, _last_variant, [&](auto &&variant) {
                return !libjst::is_insertion(variant) || libjst::position(variant) != branch_position;
            });

            // second: if next variant is not already the next valid we need to search it.
            if (_next_variant != _last_variant && branch_end > libjst::position(*_next_variant))
            {
                _next_variant = std::ranges::lower_bound(_next_variant, _last_variant, branch_end, std::less<>{},
                                                         [] (auto &&variant) { return libjst::position(variant); });
            }

            if (_next_variant != _last_variant) {
                _state.set_range(branch_position + _offset,
                                 branch_position + _offset + std::ranges::ssize(libjst::insertion(*current_variant)) +
                                 libjst::position(*_next_variant) - branch_end);
                // _next = libjst::position(*current_variant) + std::ranges::size(libjst::insertion(*current_variant)) +
                //                 libjst::position(*_next_variant) - branch_end;
            }
            _offset += std::ranges::ssize(libjst::insertion(*current_variant)) - libjst::deletion(*current_variant);
        }

        value_type const & operator*() const noexcept
        {
            return _state;
        }

        bool has_next() const noexcept
        {
            return _next_variant != _last_variant && (libjst::position(*_next_variant) + _offset) < _max_branch_size;
        }

        jst_node_branch next()
        {
            assert(_next_variant != _last_variant);
            jst_node_branch branch_node{_state, _next_variant, _last_variant, _max_branch_size, _offset};

            auto const branch_position = libjst::position(_next_variant) + _offset; // the position of this branch
            _state.unset(libjst::coverage(*_next_variant));

            // now here we don't know, if there is something better?
            auto next_position = (++_next_variant != _last_variant)
                               ? libjst::position(_next_variant) + _offset
                               : _max_branch_size;
            assert(next_position >= branch_position);
            _state.set_range(branch_position, next_position);

            // if (_value.coverage().any()) { // is there at least one sequence covering the base branch at this site.
            //     if (++_next_variant != _last_variant) { // might end before.
            //         _next += libjst::position(*_next_variant) - libjst::position(prev_variant);
            //     } else {
            //         _next = _last; // set to last.
            //     }
            // } else {
            //     _next = _last; // set to last, so it is over! empty?
            // }
        }

        constexpr operator bool() const noexcept
        {
            return has_value();
        }
    };

} // namespace libjst
