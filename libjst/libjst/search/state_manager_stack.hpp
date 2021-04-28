// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::search_state_manager_stack.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <stack>

#include <libjst/search/state_manager_concept.hpp>

namespace libjst
{

/*!\brief A search state manager handling states that are popped/pushed on a stack.
 * \implements libjst::search_state_manager
 * \implements libjst::search_stack_observer
 *
 * \tparam search_state_t The type of the search state; must model std::semiregular.
 *
 * \details
 *
 * This manager is useful for searches over the libjst::journaled_sequence_tree which implements a stack based
 * traversal. The stateful search algorithms can thus keep track of the current traversal by pushing/popping their
 * state onto a stack whenever the traversal notifies about such an event.
 */
template <std::semiregular search_state_t>
class search_state_manager_stack
{
private:
    //!\brief The inner state stack.
    std::stack<search_state_t, std::vector<search_state_t>> _state_stack{{search_state_t{}}};
public:
    //!\brief The type of the managed search states.
    using state_type = search_state_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr search_state_manager_stack() = default; //!< Default.
    constexpr search_state_manager_stack(search_state_manager_stack const &) = default; //!< Default.
    constexpr search_state_manager_stack(search_state_manager_stack &&) = default; //!< Default.
    constexpr search_state_manager_stack & operator=(search_state_manager_stack const &) = default; //!< Default.
    constexpr search_state_manager_stack & operator=(search_state_manager_stack &&) = default; //!< Default.
    ~search_state_manager_stack() = default; //!< Default.
    //!\}

    //!\copydoc libjst::search_state_manager::state()
    constexpr state_type & state() noexcept
    {
        assert(!_state_stack.empty());
        return _state_stack.top();
    }

    //!\copydoc libjst::search_state_manager::state() const
    constexpr state_type const & state() const noexcept
    {
        assert(!_state_stack.empty());
        return _state_stack.top();
    }

    //!\copydoc libjst::search_stack_observer::on_push()
    constexpr void on_push()
    {
        _state_stack.push(state());
    }

    //!\copydoc libjst::search_stack_observer::on_pop()
    constexpr void on_pop()
    {
        assert(!_state_stack.empty());
        _state_stack.pop();
    }
};
}  // namespace libjst
