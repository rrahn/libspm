// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides .
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <stack>

#include <libjst/traversal/stack_publisher.hpp>
namespace libjst
{
    // Can we add subscription to the tree directly!?
    template <typename tree_t>
    class tree_traverser_base : public stack_publisher {
    private:
        using node_type = libjst::tree_node_t<tree_t>;
        using branch_type = std::stack<node_type>;

        std::reference_wrapper<tree_t const> _tree;
        branch_type _branch{};

        class sentinel;
        class iterator;

    public:
        explicit tree_traverser_base(tree_t const & tree) noexcept : _tree{std::cref(tree)}
        {}
        explicit tree_traverser_base(tree_t && tree) noexcept = delete;

        constexpr iterator begin() noexcept {
            return iterator{*this};
        }

        sentinel end() noexcept {
            return sentinel{*this};
        }
    };

    template <typename tree_t>
    class tree_traverser_base<tree_t>::iterator {
    private:

        friend tree_traverser_base;

        tree_traverser_base * _host{};

        explicit iterator(tree_traverser_base & host) : _host{std::addressof(host)}
        {
            visit_next(libjst::root(_host->_tree.get()));
        }
    public:

        using value_type = libjst::node_label_t<node_type>;
        using reference = libjst::node_label_t<node_type>;
        using difference_type = std::ptrdiff_t;
        using pointer = std::conditional_t<std::is_lvalue_reference_v<reference>, value_type *, void>;
        using iterator_category = std::input_iterator_tag;

        iterator() = default;

        constexpr reference operator*() const noexcept {
            return *active_node();
        }

        constexpr auto operator->() const noexcept -> pointer
            requires (!std::same_as<pointer, void>) {
            return std::addressof(*this);
        }

        constexpr iterator & operator++() {
            // we can do two basic operations.
            node_type & parent = active_node();
            auto alt_child = parent.next_alt();
            auto ref_child = parent.next_ref(); // can be moved and yield selfupdate.

            if (alt_child && ref_child) {
                parent = std::move(*ref_child);
                visit_next(std::move(*alt_child));
            } else if (alt_child || ref_child) {
                parent = (ref_child) ? std::move(*ref_child) : std::move(*alt_child);
            } else {
                backtrack();
            }
            return *this;
        }

    private:

        constexpr friend bool operator==(iterator const & lhs, sentinel const &) noexcept {
            return lhs.branch().empty();
        }

        constexpr node_type const & active_node() const noexcept {
            assert(!branch().empty());
            return branch().top();
        }

        constexpr node_type & active_node() noexcept {
            assert(!branch().empty());
            return branch().top();
        }

        constexpr void visit_next(node_type && new_node) noexcept {
            branch().push(std::move(new_node));
            _host->notify_push();
        }

        constexpr void backtrack() noexcept {
            assert(!branch().empty());
            branch().pop();
            _host->notify_pop();
        }

        constexpr branch_type & branch() const noexcept {
            assert(_host != nullptr);
            return _host->_branch;
        }
    };

    template <typename tree_t>
    class tree_traverser_base<tree_t>::sentinel {
    private:
        friend tree_traverser_base;

        using sink_type = libjst::tree_sink_t<tree_t const>;

        sink_type _sink{};

        constexpr sentinel(tree_traverser_base const & host) noexcept : _sink{libjst::sink(host._tree.get())}
        {}

    public:
        constexpr sentinel() = default;

    private:

        constexpr friend bool operator==(sentinel const & lhs, node_type const & rhs) {
            return lhs._sink == rhs;
        }
    };

}  // namespace libjst
