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

#include <libjst/structure/concept_jst.hpp>
#include <libjst/concept.hpp>
#include <libjst/journal.hpp>
#include <libjst/traversal/variant_branch_node.hpp>
#include <libjst/variant/variant_proxy_null.hpp>

namespace libjst
{
    template <typename label_t, typename variant_iterator_t>
    class source_branch_node
    {
    private:

        using variant_t = std::iter_value_t<variant_iterator_t>;
        using position_t = libjst::variant_position_t<variant_t>;

        // variant_store_node_t
        label_t _label{}; // the one we need to create.
        variant_iterator_t _next_variant{};
        variant_iterator_t _last_variant{};
        position_t _context_size{};

    public:
        // what kind of node type?
        using branch_node_type = variant_branch_node<label_t, variant_iterator_t>;
        using label_type = label_t;

        source_branch_node() = default;
        explicit source_branch_node(label_t label,
                                    variant_iterator_t next_variant,
                                    variant_iterator_t last_variant,
                                    position_t const context_size) :
            _label{std::move(label)},
            _next_variant{std::move(next_variant)},
            _last_variant{std::move(last_variant)},
            _context_size{context_size - 1}
        {
            assert(context_size > 0);

            _label.reset(null_variant{0, _label.coverage()}, _next_variant);


            //     _next = libjst::position(*_next_variant);
            //     _last = _next + std::ranges::size(libjst::insertion(*_next_variant)) + _context_size;
            // } else {
            //     _next = std::ranges::size(_label.sequence());
            //     _last = _next;
            // }
        }

        label_type const & operator*() const noexcept
        {
            return _label;
        }

        label_type const * operator->() const noexcept
        {
            return std::addressof(_label);
        }

        bool is_leaf() const noexcept
        {
            return _next_variant == _last_variant;
        }

        branch_node_type branch()
        {
            assert(!is_leaf());

            position_t const branch_size = std::ranges::ssize(libjst::insertion(*_next_variant)) + _context_size;
            branch_node_type variant_root{_label, _next_variant, _last_variant, branch_size};

            position_t const label_begin = libjst::position(_next_variant);
            update_label(null_variant{label_begin, _label.coverage()}, ++_next_variant);

            // move to the next
            // ++_next_variant; // move to the next variant
            // if (++_next_variant != _last_variant)
            //     _label.update(*_next_variant);  // or update, sufficient with one parameter?

            // if (_next_variant != _last_variant) {
            //     _next = libjst::position(*_next_variant);
            //     _last = _next + std::ranges::size(libjst::insertion(*_next_variant)) + _context_size;
            // } else {
            //     _next = std::ranges::size(_label.sequence());
            //     _last = std::ranges::size(_label.sequence());
            // }
            return variant_root;
        }

        // returns false if label has no value!
        constexpr operator bool() const noexcept
        {
            return _label.has_value();  // if (!node) -> node is nil
        }

    private:

        template <typename null_variant_t, typename variant_t>
        void update_label(null_variant_t && null_variant, variant_t const & variant)
        {
            position_t label_size = ((variant != _last_variant) ?
                                            libjst::position(variant) :
                                            std::numeric_limits<position_t>::max()) - libjst::position(null_variant);

            _label.reset((null_variant_t &&) null_variant, label_size);

        }
    };
} // namespace libjst
