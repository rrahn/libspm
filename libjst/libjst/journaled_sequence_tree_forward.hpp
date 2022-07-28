// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::journaled_sequence_tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <optional>
#include <ranges>
#include <stack>
#include <type_traits>

#include <seqan3/utility/detail/multi_invocable.hpp>

#include <libjst/journal.hpp>

namespace libjst {

// JST is a data structure and we may add features by wrapping it.
// but now we want to search it with some algorithm.
// this depends on the structure type.
namespace _forward {


template <typename jst_t, typename searcher_t, typename receiver_t>
struct operation {
private:
    using search_state_t = std::decay_t<searcher_t>;
    using node_t = typename jst_t::node_type;

    static_assert(std::semiregular<search_state_t>);

    using algorithm_stack_t = std::stack<search_state_t, std::vector<search_state_t>>;
    using node_stack_t = std::stack<node_t, std::vector<node_t>>;

    jst_t const & _jst;
    searcher_t _searcher;
    receiver_t _receiver;

    // algorithm to start the operation
    void start() {
        algorithm_stack_t algorithm_stack{};
        node_stack_t node_stack{};

        algorithm_stack.push(std::move(_searcher));
        node_stack.emplace(_jst, _searcher.window_size());

        while (!node_stack.empty()) {
            search_state_t & algorithm = algorithm_stack.top();
            node_t & node = node_stack.top();
            // now every time the algorithm is invoked we need
            algorithm(node.sequence(), [&] (auto && hit) { _receiver.set_next(std::forward<decltype(hit)>(hit)); });
            if (!node.at_end()) {
                auto [branch, split] = bifurcate(std::move(node));
                assert(branch.has_value() || split.has_value());

                if (split.has_value()) { // move assign the split node, keep algorithm state.
                    node_stack.top() = std::move(*split);
                    if (branch.has_value()) {
                        node_stack.push(std::move(*branch)); // move the variant node to the top of the stack.
                        algorithm_stack.push(algorithm); // copy the algorithm state.
                    }
                } else { // split node was not added, so we can reuse the memory location on the stack.
                    assert(branch.has_value());
                    node_stack.top() = std::move(*branch);
                }
            }
        }
        std::forward<receiver_t>(_receiver).set_value(); // now we remove the handle of the receiver
    }
};

template <typename jst_t, typename searcher_t>
struct sender {

    jst_t const & _jst;
    searcher_t _searcher;

    template <typename receiver_t>
    auto connect(receiver_t && receiver) {
        using operation_t = operation<jst_t, searcher_t, receiver_t>;
        return operation_t{_jst, std::forward<searcher_t>(_searcher), std::forward<receiver_t>(receiver)};
    }
};
} // namespace _forward

template <typename jst_t>
class journaled_sequence_tree_forward : protected jst_t {
private:

    using typename jst_t::position_type;
    using typename jst_t::sequence_type;
    using typename jst_t::coverage_type;
    using typename jst_t::branch_events_type;

public:
    class node_type {

        enum branch_kind {
            base,
            variant
        };

        using variant_iterator = std::ranges::iterator_t<branch_events_type const &>;
        using journal_type = journal<position_type, sequence_type>;

        journal_type _journal{};
        variant_iterator _next_variant{};
        coverage_type _coverage{};
        size_t _first{};
        size_t _next{};
        size_t _last{};
        size_t _window_size{};
        journaled_sequence_tree_forward const * _jst{};
        branch_kind _kind{branch_kind::base};

    public:

        node_type() = default;
        explicit node_type(journaled_sequence_tree_forward const & jst, size_t window_size) :
            _jst{std::addressof(jst)},
            _journal{_jst->reference()},
            _next_variant{std::ranges::begin(_jst->_branch_event_queue)},
            _coverage{coverage_type(_jst->size(), true)},
            _window_size{window_size - 1}
        {
            assert(window_size > 0);

            if (_next_variant != std::ranges::end(_jst->_branch_event_queue)) {
                auto && next_variant = (*_next_variant).event_handle();
                _next = next_variant.position();
                _last = std::min<size_t>(_next + next_variant.insertion_size() + _window_size, _jst->reference().size());
            }
        }

        auto sequence() const noexcept {
            auto seq = _journal.sequence();
            auto it = seq.begin() + _first;
            return std::ranges::subrange{it, it + (_next - _first)};
        }

        bool at_end() const noexcept {
            return _next == _last;
        }

        friend auto bifurcate(node_type && parent)
            -> std::pair<std::optional<node_type>, std::optional<node_type>> noexcept {
            // here we need to split into the different nodes.
            // this is the function that becomes important in our scenario.
            // ------------------
            // create branch node

            std::optional<node_type> branch_node{std::nullopt_t};
            node_type child{};
            child._jst = parent._jst;
            child._next_variant = parent._next_variant; // base variant
            child._coverage = parent._coverage & (*child._next_variant).coverage();
            if (child._coverage.any()) {
                child._window_size = parent._window_size;
                child._first = parent._next;
                child._last = parent._last;
                child._journal = parent._journal;
                child._kind = branch_kind::variant;
                record_sequence_variant(child, *(*child._next_variant).event_handle());
                // find first branch candidate:
                // First find first variant that is not an insertion including the current variant.
                child._next_variant = std::ranges::find_if(++child._next_variant,
                                                           std::ranges::end(child._jst->_branch_event_queue),
                                                           [] (auto && it) {
                    return !(*it).event_handle()->is_insertion();
                });
                // second: if next variant is not already the next valid we need to search it.
                auto last_position = [] (auto && branch_event) {
                    return branch_event.position() + branch_event.event_handle()->deletion_size();
                };
                if (child._next_variant != std::ranges::end(child._jst->_branch_event_queue) &&
                    last_position(*parent._next_variant) > (*child._next_variant).position()) {
                    child._next_variant =
                        std::ranges::upper_bound(child._next_variant,
                                                 std::ranges::end(child._jst->_branch_event_queue),
                                                 last_position(*child._next_variant),
                                                 [] (auto && branch_event) { return branch_event.position(); };);
                }
                // update next relative position for bifurcation.
                size_t const next_bifurcation = parent._next +
                                                (*child._next_variant).position() -
                                                last_position(*parent._next_variant) +
                                                (*child._next_variant).event_handle()->insertion_size();
                child._next = std::min(next_bifurcation, child._last);
                branch_node = std::move(child);
            }

            // ------------------
            // create split node

            std::optional<node_type> split_node{std::nullopt_t};
            parent._first = parent._next;
            auto & last_variant = *parent._next_variant;
            ++parent._next_variant;
            if (parent._kind == branch_kind::base) { // update base branch node
                parent._next = parent._jst->reference().size();
                size_t insert_size = parent._window_size;
                if (parent._next_variant != std::ranges::end(parent._jst->_branch_event_queue)) {
                    parent._next = (*parent._next_variant).position();
                    insert_size += (*parent._next_variant).event_handle()->insertion_size();
                }
                parent._last = std::min(parent._next + insert_size, parent._jst->reference().size());
                split_node = std::move(parent);
            } else { // update variant branch node
                parent._coverage = parent._coverage.and_not(last_variant.coverage());
                if (parent._coverage.any()) {
                    // can point to last one.
                    size_t next_position = parent._jst->reference().size();
                    if (parent._next_variant != std::ranges::end(parent._jst->_branch_event_queue)) {
                        next_position = (*parent._next_variant).position();
                    }
                    assert(last_variant.position() <= next_position);

                    parent._next = std::min(parent._next + last_variant.position() - next_position, parent._last);

                    split_node = std::move(parent);
                }
            }
            return {std::move(branch_node), std::move(split_node)};
        }

    private:
        template <typename sequence_variant_t>
        static void record_sequence_variant(node_type & node, sequence_variant_t const & variant) noexcept
        {
            using insertion_t = sequence_variant_t::insertion_type;
            using deletion_t = sequence_variant_t::deletion_type;

            std::visit([&] (auto & delta_kind) {
                seqan3::detail::multi_invocable {
                    [&] (insertion_t const & e) { node._journal.record_insertion(node._first, e.value()); },
                    [&] (deletion_t const & e) { node._journal.record_deletion(node._first, e.value()); },
                    [&] (auto const & e) { node._journal.record_substitution(node._first, e.value()); }
                } (delta_kind);
            }, variant.delta_variant());
        }
    };

    journaled_sequence_tree_forward() = delete;
    explicit journaled_sequence_tree_forward(jst_t && jst) noexcept : jst_t{std::move(jst)}
    {}

    // return sender!
    template <typename searcher_t>
    auto search(searcher_t && searcher) const & {
        return _forward::sender<journaled_sequence_tree_forward, searcher_t>{*this, std::forward<searcher_t>(searcher)};
    }

    // forbid rvalue invocation
    template <typename searcher_t>
    void search(searcher_t) && = delete;
};
// now we can have the algorithm/sender here.

} // namespace libjst
