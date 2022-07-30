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

enum struct resume_traversal {
    head_on_breakpoint,
    tail_on_breakpoint
};

// JST is a data structure and we may add features by wrapping it.
// but now we want to search it with some algorithm.
// this depends on the structure type.
namespace _forward {


template <typename jst_t, typename searcher_t, typename receiver_t>
struct operation {

    using search_state_t = std::decay_t<searcher_t>;
    using node_t = typename jst_t::node_type<std::remove_cvref_t<searcher_t>::resume_policy>;

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
            } else {
                algorithm_stack.pop();
                node_stack.pop();
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

    // using typename jst_t::position_type;
    using position_type = uint32_t;
    using typename jst_t::sequence_type;
    using typename jst_t::coverage_type;
    using typename jst_t::branch_events_type;

public:

    template <resume_traversal resume_policy = resume_traversal::tail_on_breakpoint>
    class node_type {

        enum branch_kind {
            base,
            variant
        };

        using variant_iterator = std::ranges::iterator_t<branch_events_type const &>;
        using journal_type = journal<position_type, sequence_type const &>;

        journal_type _journal{};
        variant_iterator _next_variant{};
        coverage_type _coverage{};
        size_t _first{};
        size_t _next{};
        size_t _last{};
        size_t _window_size{};
        branch_events_type const * _variants{};
        sequence_type const * _base_sequence{};
        branch_kind _kind{branch_kind::base};

    public:

        node_type() = default;
        explicit node_type(journaled_sequence_tree_forward const & jst, size_t window_size) :
            _journal{jst.reference_at(0)},
            _next_variant{std::ranges::begin(jst.variants())},
            _coverage{coverage_type(jst.size(), true)},
            _window_size{window_size - 1},
            _variants{std::addressof(jst.variants())},
            _base_sequence{std::addressof(jst.reference_at(0))}
        {
            assert(window_size > 0);

            if (_next_variant != std::ranges::end(variants())) {
                auto && next_variant = *(*_next_variant).event_handle();
                _next = next_variant.position().offset;
                _last = std::min<size_t>(_next + next_variant.insertion_size() + _window_size, base_sequence().size());
            } else {
                _next = base_sequence().size();
                _last = _next;
            }
        }

        auto sequence() const noexcept {
            auto seq = _journal.sequence();
            size_t head_position = _first;
            if constexpr (resume_policy == resume_traversal::tail_on_breakpoint) {
                head_position -= std::min(_window_size, _first);
            }
            auto it = seq.begin() + head_position;
            return std::ranges::subrange{it, it + (_next - head_position)}; // borrowed range
        }

        bool at_end() const noexcept {
            return _next == _last;
        }

        template <typename node_t>
            requires std::same_as<std::remove_reference_t<node_t>, node_type>
        friend auto bifurcate(node_t && parent) noexcept
            -> std::pair<std::optional<std::remove_reference_t<node_t>>, std::optional<std::remove_reference_t<node_t>>>
        {
            // here we need to split into the different nodes.
            // this is the function that becomes important in our scenario.
            // ------------------
            // create branch node

            std::optional<node_type> branch_node{std::nullopt};
            node_type child{};
            child._variants = parent._variants;
            child._base_sequence = parent._base_sequence;
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
                                                           std::ranges::end(child.variants()),
                                                           [] (auto && branch_event) {
                    return !branch_event.event_handle()->is_insertion();
                });
                // second: if next variant is not already the next valid we need to search it.
                auto last_position = [] (auto && branch_event) {
                    return branch_event.position().offset + branch_event.event_handle()->deletion_size();
                };
                if (child._next_variant != std::ranges::end(child.variants()) &&
                    last_position(*parent._next_variant) > (*child._next_variant).position().offset) {
                    child._next_variant =
                        std::ranges::lower_bound(child._next_variant,
                                                 std::ranges::end(child.variants()),
                                                 last_position(*parent._next_variant),
                                                 std::less<>{},
                                                 [] (auto && branch_event) { return branch_event.position().offset; });
                }
                // update next relative position for bifurcation.
                size_t const next_bifurcation = parent._next +
                                                child.next_variant_position() - last_position(*parent._next_variant) +
                                                (*parent._next_variant).event_handle()->insertion_size();
                child._next = std::min(next_bifurcation, child._last);
                branch_node = std::move(child);
            }

            // ------------------
            // create split node

            std::optional<node_type> split_node{std::nullopt};
            parent._first = parent._next;
            auto & last_variant = *parent._next_variant;
            ++parent._next_variant;
            if (parent._kind == branch_kind::base) { // update base branch node
                parent._next = parent.base_sequence().size();
                size_t insert_size = parent._window_size;
                if (parent._next_variant != std::ranges::end(child.variants())) {
                    parent._next = (*parent._next_variant).position().offset;
                    insert_size += (*parent._next_variant).event_handle()->insertion_size();
                }
                parent._last = std::min(parent._next + insert_size, parent.base_sequence().size());
                split_node = std::move(parent);
            } else { // update variant branch node
                parent._coverage = parent._coverage.and_not(last_variant.coverage());
                if (parent._coverage.any()) {
                    // can point to last one.
                    size_t next_position = parent.base_sequence().size();
                    if (parent._next_variant != std::ranges::end(child.variants())) {
                        next_position = (*parent._next_variant).position().offset;
                    }
                    assert(last_variant.position().offset <= next_position);

                    parent._next = std::min(parent._next + next_position - last_variant.position().offset, parent._last);

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

            // TODO: Revert dirty fix and let journaled sequence tree decide which type of sequence to allocate.
            std::visit([&] (auto & delta_kind) {
                seqan3::detail::multi_invocable {
                    [&] (insertion_t const & e) { node._journal.record_insertion(node._first, std::string_view{e.value().data(), e.value().size()}); },
                    [&] (deletion_t const & e) { node._journal.record_deletion(node._first, e.value()); },
                    [&] (auto const & e) { node._journal.record_substitution(node._first, std::string_view{e.value().data(), e.value().size()}); }
                } (delta_kind);
            }, variant.delta_variant());
        }

        branch_events_type const & variants() const noexcept {
            return *_variants;
        }

        sequence_type const & base_sequence() const noexcept {
            return *_base_sequence;
        }

        size_t next_variant_position() const noexcept {
            return (_next_variant == std::ranges::end(variants())) ? base_sequence().size()
                                                                   : (*_next_variant).position().offset;
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

protected:

    branch_events_type const & variants() const noexcept {
        return this->_branch_event_queue;
    }
};
// now we can have the algorithm/sender here.

} // namespace libjst
