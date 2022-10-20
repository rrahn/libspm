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

#include <iterator>
#include <ranges>

#include <libjst/structure/concept_jst.hpp>
#include <libjst/variant/variant_proxy_offset.hpp>
#include <libjst/variant/variant_proxy_coverage_transform.hpp>

namespace libjst
{
    // Clearly needs some refactoring!
    template <typename label_t, typename variant_iterator_t>
    class variant_branch_node
    {
    private:
        using position_t = variant_position_t<std::iter_value_t<variant_iterator_t>>;

        label_t _label{};
        variant_iterator_t _next_variant{};
        variant_iterator_t _last_variant{};
        position_t _remaining_size{};
        position_t _offset{};

    public:

        using label_type = label_t;

        variant_branch_node() = default;
        variant_branch_node(variant_branch_node const &) = default;
        variant_branch_node(variant_branch_node &&) = default;
        variant_branch_node & operator=(variant_branch_node const &) = default;
        variant_branch_node & operator=(variant_branch_node &&) = default;

        // construction from base node
        explicit variant_branch_node(label_t parent_label, // parent_label
                                 variant_iterator_t const current_variant,
                                 variant_iterator_t last_variant,
                                 position_t const remaining_size,
                                 position_t const offset = 0) :
            _label{std::move(parent_label)},
            _last_variant{std::move(last_variant)},
            _remaining_size{remaining_size},
            _offset{offset}
        {
            assert(_next_variant != _last_variant);

            _next_variant = find_next(current_variant);
            update_label(current_variant, _next_variant, [&] (coverage_t const & branch_coverage) {
                return std::move(_label.coverage()) & branch_coverage;
            });

            // position_t next_position = ((_next_variant != _last_variant) ?
            //                                 libjst::position(*_next_variant) :
            //                                 std::numeric_limits<position_t>::max());
            // position_t label_size = min(next_position - libjst::position(*current_variant), _remaining_size);

            // _label.reset(coverage_transform_variant{offset_variant{*current_variant, offset},
            // }, label_size);
            // _remaining_size -= label_size;
            _offset += std::ranges::ssize(libjst::insertion(*current_variant)) - libjst::deletion(*current_variant);
        }

        label_type const & operator*() const noexcept
        {
            return _label;
        }

        bool is_leaf() const noexcept
        {
            return _next_variant != _last_variant || _remaining_size == 0;
        }

        variant_branch_node branch()
        {
            assert(_next_variant != _last_variant);
            variant_iterator_t current_variant = _next_variant;
            variant_branch_node branch_node{_label, current_variant, _last_variant, _remaining_size, _offset};

            update_label(current_variant, ++_next_variant, [&] (coverage_t const & branch_coverage) -> coverage_t {
                coverage_t tmp = std::move(_label.coverage());
                return tmp.and_not(branch_coverage);
            });

            // position_t next_position = ((++_next_variant != _last_variant) ?
            //                                 libjst::position(*_next_variant) :
            //                                 std::numeric_limits<position_t>::max());
            // position_t label_size = min(next_position - libjst::position(*current), _remaining_size);
            // _label.reset(coverage_transform_variant{offset_variant{*current, _offset},
            // }, label_size);
            // _remaining_size -= label_size;

//   - dist_to_next = (++next != end) : pos(next) - pos(curr) ? rest
//   - label_size = min(dist_to_next, rest)
//   - s.reset(invert_coverage{*next, off}, label_size) // nullifies the coverage
//   - rest -= label_size

            // auto const branch_position = libjst::position(_next_variant) + _offset; // the position of this branch
            // _label.unset(libjst::coverage(*_next_variant));

            // // now here we don't know, if there is something better?
            // auto next_position = (++_next_variant != _last_variant)
            //                    ? libjst::position(_next_variant) + _offset
            //                    : _remaining_size;
            // assert(next_position >= branch_position);
            // _label.set_range(branch_position, next_position);



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
            return _label.has_value();  // if (!node) -> node is nil
        }

        private:

        template <typename fn_t>
        void update_label(variant_iterator_t current_branch,
                          variant_iterator_t next_branch,
                          fn_t && coverage_fn)
        {
            position_t next_position = ((next_branch != _last_variant) ?
                                            libjst::position(*next_branch) :
                                            std::numeric_limits<position_t>::max());
            position_t const label_size = min(next_position - libjst::position(*current_branch), _remaining_size);

            // but we need to wrap something here as well.
            null_variant nullvar{libjst::position(*current_branch), libjst::coverage(*current_branch)};
            _label.reset(coverage_transform_variant{offset_variant{std::move(nullvar), _offset}, (fn_t &&) coverage_fn},
                         label_size);
            _remaining_size -= label_size;
        }

        variant_iterator find_next(variant_iterator it) const
        {
            position_t const branch_position = libjst::position(*it);
            position_t const branch_end = branch_position + libjst::deletion(*it);
            it = std::ranges::find_if(++it, _last_variant, [&] (auto &&variant) {
                return !libjst::is_insertion(variant) || libjst::position(variant) != branch_position;
            });

            // second: if next variant is not already the next valid we need to search it.
            if (it != _last_variant && branch_end > libjst::position(*it)) {
                it = std::ranges::lower_bound(it, _last_variant, branch_end, std::less<>{},
                                                         [] (auto &&variant) { return libjst::position(variant); });
            }

            return it;
        }
    };

} // namespace libjst
