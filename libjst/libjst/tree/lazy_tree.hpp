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

#include <libjst/set/concept_set.hpp>
#include <libjst/tree/stack_publisher.hpp>

namespace libjst
{
    // now what are we expecting to do here?
    // and who implements the actual search?
    // template <typename jst_t, typename searcher_t, typename receiver_t>
    template <typename node_t> // the node is customisable -> requires std::semiregular and ?
    class lazy_tree : protected stack_publisher// offer range interface as stream, something special?
    {
    private:
        using node_stack_t = std::stack<node_t, std::vector<node_t>>;
        class iterator;
        class sentinel;

        node_stack_t _current_branch{};

    public:

        size_t _prune_count{};
        size_t _branch_count{};

        template <observable_stack ...subscriber_t>
        explicit lazy_tree(subscriber_t &...subscriber) noexcept
        {
            (stack_publisher::subscribe(subscriber), ...);
        }

        template <typename _node_t, observable_stack ...subscriber_t>
            requires (std::same_as<std::remove_cvref_t<_node_t>, node_t>)
        explicit lazy_tree(_node_t &&root, subscriber_t &...subscriber) noexcept :
            _current_branch{(_node_t &&)root}
        {
            (stack_publisher::subscribe(subscriber), ...);
            stack_publisher::notify_push();
        }

        // template <traversable_journaled_sequence_tree jst_t, observable_stack ...subscriber_t>
        // explicit lazy_tree(jst_t const & jst, size_t const window_size, subscriber_t &...subscriber) noexcept :
        //     _window_size{window_size}
        // {
        //     // requires jst_node to be instantiatable with jst itself and a window size.
        //     (stack_publisher::subscribe(subscriber), ...);
        //     _current_branch.emplace(jst, _window_size);
        //     stack_publisher::notify_push();
        // }

        void reset(node_t node) {
            // normally we need to reset everything also, the subscription!
            assert(_current_branch.empty());
            _current_branch.push(std::move(node));
            stack_publisher::notify_push();
        }

        // input iterator
        iterator begin() noexcept
        {
            return iterator{this};
        }

        iterator begin() const noexcept = delete;

        sentinel end() noexcept
        {
            return {};
        }

        sentinel end() const noexcept = delete;
        // algorithm to start the operation
    };

    template <typename node_t, typename ...subscriber_t>
    lazy_tree(node_t &&, subscriber_t &&...) -> lazy_tree<std::remove_cvref_t<node_t>>;

    template <typename node_t>
    class lazy_tree<node_t>::iterator
    {
    private:
        friend sentinel;
        lazy_tree * _tree{};
    public:

        using value_type = node_t;
        using reference = node_t const &;
        using pointer = node_t const *;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::input_iterator_tag;

        iterator() = default;
        iterator(lazy_tree * tree) : _tree{tree}
        {}
        iterator(iterator const &) = delete;
        iterator(iterator &&) = default;
        iterator &operator=(iterator const &) = delete;
        iterator &operator=(iterator &&) = default;

        reference operator*() const noexcept
        {
            assert(_tree != nullptr);
            return current_node();
        }

        pointer operator->() const noexcept
        {
            assert(_tree != nullptr);
            return std::addressof(current_node());
        }

        iterator & operator++() noexcept
        {
            assert(_tree != nullptr);
            node_t & parent = current_node();
            if (parent.has_value())
            {
                auto child = parent.split();
                // auto [branch, split] = bifurcate(std::move(node)); // remove the top node element
                assert(child.has_value() || parent.has_value()); // optional values
                _tree->_prune_count += !child.has_value() + !parent.has_value();
                _tree->_branch_count += child.has_value() + parent.has_value();
                if (parent.has_value()) // if split has no value then branch state must have?
                {
                    if (child.has_value()) {
                        _tree->_current_branch.push(std::move(*child));
                        _tree->notify_push();
                    }
                }
                else
                { // split node was not added, so we can reuse the memory location on the stack.
                    assert(child.has_value());
                    current_node() = std::move(*child);
                }
            }
            else
            {
                _tree->_current_branch.pop(); // sets to the next node!
                _tree->notify_pop();
            }
            return *this;
        }

        node_t operator++(int) noexcept(std::is_nothrow_copy_constructible_v<node_t>)
        {
            assert(_tree != nullptr);
            node_t tmp{*(*this)};
            ++(*this);
            return tmp;
        }
    private:
        node_t & current_node() noexcept
        {
            return _tree->_current_branch.top();
        }

        node_t const & current_node() const noexcept
        {
            return _tree->_current_branch.top();
        }
    };

    template <typename node_t>
    class lazy_tree<node_t>::sentinel
    {
    public:
        constexpr bool operator==(iterator const & rhs) const noexcept
        {
            return rhs._tree == nullptr || rhs._tree->_current_branch.empty();
        }
    };

    // input is searcher and jst
    // node type depends on the tree type.
    // so depending on the search and tree we can get different implementations.
    // simply check if the searcher is stateful!
    // lazy_tree<resumable_node> jst_tree{jst, window_size(searcher)}; // other tree type
    // state_tree state_tree{libjst::state(searcher)}; // if stateful
    // linked_tree tree{std::move(jst_tree), std::move(state_tree)}; // composite!
    // for (auto && [node, searcher] : tree) // when to drop the algorithm? // => while (!node_stack.empty()) // always one step
    // {
    //     searcher(node.label(), [] (auto & hit) { _reveiver.set_next(wrap_hit(hit)); }  ); // so we can wrap the searcher/ we can wrap the receiver to customise the information.
    // }

    // how can we implement the search as part of some traversal algorithm?

    //  void start()
    //     {
    //         algorithm_stack_t algorithm_stack{}; // this is handled somewhere else!
    //         lazy_tree tree{window_size(searcher)};
    //         node_stack_t node_stack{};

    //         search_t op = libjst::search_operation((searcher_t &&) _searcher);
    //         node_stack.emplace(_jst, libjst::window_size(op));
    //         algorithm_stack.push(std::move(op));

    //         while (!node_stack.empty())
    //         {
    //             search_t &algorithm = algorithm_stack.top();
    //             node_t &node = node_stack.top();
    //             // now every time the algorithm is invoked we need
    //             algorithm(node.sequence(), [&](auto &&hit)
    //                       { _receiver.set_next(std::forward<decltype(hit)>(hit)); });
    //             if (!node.at_end())
    //             {
    //                 auto [branch, split] = bifurcate(std::move(node));
    //                 assert(branch.has_value() || split.has_value());

    //                 if (split.has_value())
    //                 { // move assign the split node, keep algorithm state.
    //                     node_stack.top() = std::move(*split);
    //                     if (branch.has_value())
    //                     {
    //                         node_stack.push(std::move(*branch)); // move the variant node to the top of the stack.
    //                         algorithm_stack.push(algorithm);     // copy the algorithm state.
    //                     }
    //                 }
    //                 else
    //                 { // split node was not added, so we can reuse the memory location on the stack.
    //                     assert(branch.has_value());
    //                     node_stack.top() = std::move(*branch);
    //                 }
    //             }
    //             else
    //             {
    //                 algorithm_stack.pop();
    //                 node_stack.pop();
    //             }
    //         }
    //         std::forward<receiver_t>(_receiver).set_value(); // now we remove the handle of the receiver
    //     }
}  // namespace libjst
