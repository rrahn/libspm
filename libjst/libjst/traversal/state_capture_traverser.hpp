// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides jst traversal by reinstantiating the state of the pattern from a previous node.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <stack>
#include <bitset>

#include <seqan3/core/debug_stream.hpp>

#include <libjst/sequence_tree/concept.hpp>
#include <libjst/sequence_tree/left_extend_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/prune_unsupported.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

namespace libjst
{

    struct state_capture_traverser {

        template <typename state_t>
        class state_manager;

        template <typename tree_t, typename pattern_t, typename callback_t>
        constexpr void operator()(tree_t && tree, pattern_t && pattern, callback_t && callback) const {
            if (libjst::window_size(pattern) == 0)
                return;

            auto search_tree = tree | libjst::labelled()
                                    | libjst::coloured()
                                    | trim(libjst::window_size(pattern) - 1)
                                    | prune_unsupported()
                                    | merge(); // make big nodes

            state_manager<pattern_t> listening_pattern{(pattern_t &&) pattern};
            tree_traverser_base traversal_path{search_tree};
            traversal_path.subscribe(listening_pattern);
            // we need to add another stack but extern of the algorithm.
            for (auto it = traversal_path.begin(); it != traversal_path.end(); ++it) {
                auto && label = *it;
                listening_pattern(label.sequence(), [&] (auto && label_it) {
                    callback(std::move(label_it), label); // either cargo offers access to node context or not!
                });
            }
        }
    };

    template <typename matcher_t>
    class state_capture_traverser::state_manager {
    private:

        using state_t = libjst::matcher_state_t<matcher_t>;
        using state_stack_t = std::stack<state_t>;

        matcher_t _matcher;
        state_stack_t _states{};

    public:

        constexpr explicit state_manager(matcher_t matcher) noexcept :
            _matcher{(matcher_t &&) matcher},
            _states{}
        {}

        template <typename ...args_t>
        constexpr auto operator()(args_t && ...args) noexcept(std::is_nothrow_invocable_v<matcher_t, args_t...>)
            -> std::invoke_result_t<matcher_t, args_t...> {
            return _matcher((args_t &&) args...);
        }

        constexpr void notify_push() {
            _states.push(_matcher.capture());
        }

        constexpr void notify_pop() {
            assert(!_states.empty());
            _matcher.restore(_states.top());
            _states.pop();
        }

        constexpr void print_state() const noexcept {
            seqan3::debug_stream << "[";
            std::ranges::for_each(_matcher.capture(), [] (uint64_t word) {
                seqan3::debug_stream << std::bitset<64>{word} << ", ";
            });
            seqan3::debug_stream << "]\n";
        }
    };

    // template <typename matcher_t>
    // state_manager(matcher_t &&) -> state_capture_traverser::state_manager<matcher_t>;

}  // namespace libjst
