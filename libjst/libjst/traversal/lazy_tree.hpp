// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a lazy tree implementation that builds the node on the fly.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <stack>
#include <vector>

#include <libjst/structure/concept_jst.hpp>
#include <libjst/traversal/stack_publisher.hpp>

namespace libjst
{
    // now what are we expecting to do here?
    // and who implements the actual search?
    // template <typename jst_t, typename searcher_t, typename receiver_t>
    template <traversable_journaled_sequence_tree jst_t> // the node is customisable -> requires std::semiregular and ?
    class lazy_tree : protected stack_publisher// offer range interface as stream, something special?
    {
    private:

        using store_t = variant_store_t<jst_t>;
        using variant_t = std::ranges::range_value_t<store_t>;
        using coverage_t = variant_coverage_t<variant_t>;
        using position_t = variant_position_t<variant_t>;
        using path_t = jst_path_t<jst_t>;
        using jst_node_t = std::ranges::range_value_t<path_t>;

        class node_type;
        class iterator;
        using sentinel = std::default_sentinel_t;

        using variant_path_stack_t = std::stack<node_type, std::vector<node_type>>;

        variant_path_stack_t _current_variant_branch{};
        node_type _current_source_node{}; // tree over the paths not the
        position_t _context_size{};

    public:

        size_t _prune_count{};
        size_t _branch_count{};

        template <traversable_journaled_sequence_tree _jst_t, observable_stack ...subscriber_t>
            requires (!std::same_as<std::remove_cvref_t<_jst_t>, lazy_tree> &&
                      std::constructible_from<jst_t, _jst_t>)
        explicit lazy_tree(_jst_t &&jst, position_t const context_size, subscriber_t &...subscriber) noexcept :
            _current_source_node{coverage_t(true, libjst::size(jst)),
                                 *std::ranges::begin(libjst::path(jst)),
                                 std::ranges::size(libjst::base_sequence(jst))},
            _context_size{context_size}
        {
            (stack_publisher::subscribe(subscriber), ...);
            stack_publisher::notify_push();
        }

        // input iterator
        iterator begin() noexcept
        {
            return iterator{this};
        }

        sentinel end() const noexcept
        {
            return {};
        }
    };

    template <traversable_journaled_sequence_tree jst_t, typename ...subscriber_t>
    lazy_tree(jst_t &&, subscriber_t &&...) -> lazy_tree<std::remove_cvref_t<jst_t>>;

    template <traversable_journaled_sequence_tree jst_t>
    class lazy_tree<jst_t>::iterator
    {
    private:
        using path_iterator_t = std::ranges::iterator_t<path_t>;

        lazy_tree * _tree{};
    public:

        using value_type = node_type;
        using reference = node_type const &;
        using pointer = node_type const *;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::input_iterator_tag;

        iterator() = default;
        iterator(lazy_tree * tree) : _tree{tree}
        {
        }

        iterator(iterator const &) = delete;
        iterator(iterator &&) = default;
        iterator &operator=(iterator const &) = delete;
        iterator &operator=(iterator &&) = default;

        reference operator*() const noexcept
        {
            assert(_tree != nullptr);
            assert(current_node().has_value());
            return current_node();
        }

        pointer operator->() const noexcept
        {
            assert(_tree != nullptr);
            assert(current_node().has_value());
            return std::addressof(current_node());
        }

        iterator & operator++() noexcept
        {
            assert(_tree != nullptr);

            // first we need to go into the branch
            if (!_tree->_current_variant_branch.empty()) {
                next_variant_branch();
            } else if (!current_node().is_leaf()) { // go into the next variant branch.
                auto & base_node = _tree->_current_source_node;
                auto it = base_node._jst_node.next_variant();
                _tree->_current_variant_branch.emplace(libjst::coverage(*it),
                                                       base_node._jst_node.alt(),
                                                       _tree->_context_size + std::ranges::ssize(libjst::insertion(*it)) - 1);
                base_node.next();
                _tree->notify_push();
            } // otherwise reached end of tree.
            return *this;
        }

        value_type operator++(int) noexcept(std::is_nothrow_copy_constructible_v<value_type>)
        {
            assert(_tree != nullptr);
            value_type tmp{current_node()};
            ++(*this);
            return tmp;
        }

        constexpr bool operator==(sentinel) const noexcept
        {
            return _tree == nullptr || is_nil();
        }
    private:

        void next_variant_branch()
        {
            node_type & parent = _tree->_current_variant_branch.top();
            if (!parent.is_leaf()) // current node
            {
                auto child = fork_child(parent);
                // update the parent node.
                parent._coverage.and_not(libjst::coverage(*parent._jst_node.next_variant()));
                parent.next();
                update_branch(parent, std::move(child));
            } else { // top node is a leaf and was processed last time.
                _tree->_current_variant_branch.pop(); // sets to the next node!
                _tree->notify_pop();
            }
        }

        void update_branch(node_type &parent, node_type child) noexcept
        {
            assert(child.has_value() || parent.has_value()); // one of them must not be nil
            _tree->_prune_count += !child.has_value() + !parent.has_value();
            _tree->_branch_count += child.has_value() + parent.has_value();

            if (parent.has_value()) { // if parent is not nil after update we keep it on the stack.
                if (child.has_value()) {  // if child is also not nil we push it on top of the stack.
                    _tree->_current_variant_branch.push(std::move(child));
                    _tree->notify_push();
                }
            } else { // if parent is nil, child can not and repalces parent node
                assert(child.has_value());
                parent = std::move(child);
            }
        }

        constexpr node_type fork_child(node_type const & parent) const
        {
            assert(!_tree->_current_variant_branch.empty());
            auto it = parent._jst_node.next_variant(); // but now a position without offset?
            return node_type{parent.coverage() & libjst::coverage(*it),
                             parent._jst_node.alt(),
                             std::max(0, static_cast<int32_t>(parent._remaining_size) -
                                         static_cast<int32_t>(parent.size()))};
        }

        constexpr bool is_nil() const noexcept
        {
            assert(_tree != nullptr);
            return _tree->_current_variant_branch.empty() && _tree->_current_source_node.is_leaf();
        }

        reference current_node() const noexcept
        {
            assert(_tree != nullptr);
            if (!_tree->_current_variant_branch.empty())
                return _tree->_current_variant_branch.top();
            else
                return _tree->_current_source_node;
        }
    };

    template <traversable_journaled_sequence_tree jst_t>
    class lazy_tree<jst_t>::node_type
    {
    private:

        friend iterator;

        using variant_t = std::ranges::range_value_t<store_t>;
        using position_t = variant_position_t<variant_t>;
        using coverage_t = variant_coverage_t<variant_t>;

        coverage_t _coverage{};
        jst_node_t _jst_node{};
        position_t _remaining_size{};

    public:

        node_type() = default;
        template <std::integral _position_t>
        node_type(coverage_t coverage, jst_node_t jst_node, _position_t const remaining_size) :
            _coverage{std::move(coverage)},
            _jst_node{std::move(jst_node)},
            _remaining_size{static_cast<position_t>(remaining_size)}
        {
        }

        auto sequence() const noexcept
        {
            return _jst_node.label();
        }

        coverage_t const & coverage() const noexcept
        {
            return _coverage;
        }

        bool has_value() const noexcept
        {
            return _coverage.any();
        }

        constexpr bool is_leaf() const noexcept
        {
            return _remaining_size == 0 || _jst_node.is_leaf(); //_current == std::ranges::end(_path);
        }

        constexpr ptrdiff_t size() const noexcept
        {
            return std::min(_jst_node.label_size(), _remaining_size);
        }

        constexpr position_t next_position() const noexcept
        {
            _jst_node.end_position();
            // return (_current != std::ranges::end(_path)) ?
            //     libjst::position(*_current) :
            //     std::ranges::ssize(_path.sequence());
        }

        constexpr operator bool() const noexcept
        {
            return has_value();
        }

    private:

        constexpr void next() noexcept
        {
            _jst_node.next();
        }
    };

}  // namespace libjst
