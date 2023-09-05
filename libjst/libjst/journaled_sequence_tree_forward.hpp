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

#include <cassert>
#include <concepts>
#include <optional>
#include <ranges>
#include <stack>
#include <type_traits>
#include <variant>

#include <seqan3/utility/detail/multi_invocable.hpp>

#include <libjst/concept.hpp>
#include <libjst/journal.hpp>

namespace libjst {

// JST is a data structure and we may add features by wrapping it.
// but now we want to search it with some algorithm.
// this depends on the structure type.
namespace _forward {


template <typename jst_t, typename searcher_t, typename receiver_t>
struct operation {

    using search_t = search_operation_t_old<searcher_t>;

    // TODO: needs to be part of the is resumable concept.
    // Query properties!
    static constexpr bool _is_resumable = libjst::is_resumable_v<search_t>;

    using node_t = typename jst_t::node_type<_is_resumable>;

    static_assert(std::semiregular<search_t>);

    // do not need to be global!
    using algorithm_stack_t = std::stack<search_t, std::vector<search_t>>;
    using node_stack_t = std::stack<node_t, std::vector<node_t>>;

    jst_t const & _jst;
    searcher_t _searcher;
    receiver_t _receiver;

    // algorithm to start the operation
    void start() {
        algorithm_stack_t algorithm_stack{};
        node_stack_t node_stack{};

        search_t op = libjst::search_operation_old((searcher_t &&)_searcher);
        node_stack.emplace(_jst, libjst::window_size(op));
        algorithm_stack.push(std::move(op));

        while (!node_stack.empty()) {
            search_t & algorithm = algorithm_stack.top();
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
    auto connect(receiver_t && receiver) /*conditional noexcept + return type definition*/
    {
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

    template <bool is_resumable = false>
    class node_type {

        enum branch_kind {
            base,
            variant
        };

        using variant_iterator = std::ranges::iterator_t<branch_events_type const &>;
        using journal_type = journal<position_type, sequence_type const &>;

        journal_type _journal{};
        variant_iterator _next_variant{};
        variant_iterator _last_variant{};
        coverage_type _coverage{};
        size_t _first{};
        size_t _next{};
        size_t _last{};
        size_t _window_size{};
        size_t _base_size{};
        branch_kind _kind{branch_kind::base};

    public:

        node_type() = default;
        explicit node_type(journaled_sequence_tree_forward const & jst, size_t window_size) :
            _journal{jst.reference_at(0)},
            _next_variant{std::ranges::begin(jst.variants())},
            _last_variant{std::ranges::end(jst.variants())},
            _coverage{coverage_type(jst.size(), true)},
            _window_size{window_size - 1},
            _base_size{std::ranges::size(jst.reference_at(0))}
        {
            assert(window_size > 0);

            if (_next_variant != _last_variant) {
                auto && next_variant = *(*_next_variant).event_handle();
                _next = next_variant.position().offset;
                _last = _next + next_variant.insertion_size() + _window_size;
                    //  _base_size + next_variant.insertion_size() - next_variant.deletion_size());
            } else {
                _next = _base_size;
                _last = _next;
            }
        }

        auto sequence() const noexcept {
            auto seq = _journal.sequence();
            // std::cout << "journal sequence: " << std::string{seq.begin(), seq.end()} << "\n";
            size_t head_position{};
            auto it = seq.begin();
            if constexpr (!is_resumable) {
                head_position = _first - std::min(_window_size, _first);
                assert(seq.size() >= head_position);
                assert(_next >= head_position);
                it += head_position;
            }
            size_t subrange_size = std::min(std::min(_next, _last), seq.size()) - head_position;
            return std::ranges::subrange{it, it + subrange_size, subrange_size}; // borrowed range
        }

        bool at_end() const noexcept {
            return _next >= _last;
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

            auto child_coverage = parent._coverage & (*parent._next_variant).coverage();
            if (child_coverage.any()) {
                node_type child{};
                child._next_variant = parent._next_variant;
                child._last_variant = parent._last_variant;
                child._coverage = std::move(child_coverage);
                child._base_size = parent._base_size;
                child._window_size = parent._window_size;
                child._first = parent._next;
                child._next = parent._last;
                child._last = parent._last;
                child._journal = parent._journal;
                child._kind = branch_kind::variant;
                record_sequence_variant(child, *(*child._next_variant).event_handle());
                // find first branch candidate:
                // First find first variant that is not an insertion including the current variant.
                child._next_variant = std::ranges::find_if(++child._next_variant,
                                                           child._last_variant,
                                                           [pivot = (*parent._next_variant).position()] (auto && branch_event) {
                    return !branch_event.event_handle()->is_insertion() || branch_event.position() != pivot;
                });
                // second: if next variant is not already the next valid we need to search it.
                auto last_position = [] (auto && branch_event) {
                    return branch_event.position().offset + branch_event.event_handle()->deletion_size();
                };
                // Of course! We have to check again, if the variant is at end.
                if (child._next_variant != child._last_variant &&
                    last_position(*parent._next_variant) > (*child._next_variant).position().offset) {
                        child._next_variant =
                            std::ranges::lower_bound(child._next_variant,
                                                     child._last_variant,
                                                     last_position(*parent._next_variant),
                                                     std::less<>{},
                                                     [] (auto && variant) { return variant.position().offset; });
                }

                if (child._next_variant != child._last_variant) {
                    child._next = parent._next + (*parent._next_variant).event_handle()->insertion_size() +
                                  (*child._next_variant).position().offset - last_position(*parent._next_variant);
                }
                branch_node = std::move(child);
            }

            // ------------------
            // create split node

            std::optional<node_type> split_node{std::nullopt};
            parent._first = parent._next;
            auto const & prev_variant = *parent._next_variant;
            ++parent._next_variant;
            if (parent._kind == branch_kind::base) { // update base branch node
                parent._next = parent._base_size;
                parent._last = parent._base_size;
                if (parent._next_variant != parent._last_variant) {
                    parent._next = (*parent._next_variant).position().offset;
                    parent._last = parent._next + (*parent._next_variant).event_handle()->insertion_size() +
                                   parent._window_size;
                                //    , parent._base_size +
                                //                 (*parent._next_variant).event_handle()->insertion_size() -
                                //                 (*parent._next_variant).event_handle()->deletion_size());
                }
                split_node = std::move(parent);
            } else { // update variant branch node
                parent._coverage = parent._coverage.and_not(prev_variant.coverage());
                if (parent._coverage.any()) { // is there at least one sequence covering the base branch at this site.
                    if (parent._next_variant != parent._last_variant) { // might end before.
                        parent._next += (*parent._next_variant).position().offset - prev_variant.position().offset;
                    } else {
                        parent._next = parent._last; // set to last.
                    }
                    // assert(prev_variant.position().offset <= next_position);
                    // parent._next = std::min(parent._next + next_position - prev_variant.position().offset, parent._last);
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
                    [&] (insertion_t const & e) { node._journal.record_insertion(node._first, e.value()/*std::string_view{e.value().data(), e.value().size()}*/); },
                    [&] (deletion_t const & e) { node._journal.record_deletion(node._first, e.value()); },
                    [&] (auto const & e) { node._journal.record_substitution(node._first, e.value()/*std::string_view{e.value().data(), e.value().size()}*/); }
                } (delta_kind);
            }, variant.delta_variant());
        }

        size_t next_variant_position() const noexcept {
            return (_next_variant == _last_variant) ? std::max(_last, _base_size) : (*_next_variant).position().offset;
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

    using jst_t::reference_at;
    using jst_t::reference;
    using jst_t::total_symbol_count;

protected:

    branch_events_type const & variants() const noexcept {
        return this->_branch_event_queue;
    }
};
// now we can have the algorithm/sender here.

} // namespace libjst
