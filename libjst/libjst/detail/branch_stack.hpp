// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::detail::branch_stack.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <vector>

namespace libjst::detail
{

/*!\brief An augmented stack to store the branches during the traversal of the libjst::journal_sequence_tree_enumerator.
 * \tparam branch_t The branch type to use; must model std::semiregular.
 * \tparam container_t The container type to store the branches (defaults to std::vector); must model
 *                     std::ranges::random_access_range and the
 *                     [sequence container concept](https://en.cppreference.com/w/cpp/named_req/SequenceContainer)
 *
 * \details
 *
 * This stack adaptor offers a regular stack interface but also some additional interfaces, which are helpful
 * for the traversal over the libjst::journal_sequence_tree, for example accessing branches at any level or
 * accessing the base branch.
 */
template <std::semiregular branch_t, std::ranges::random_access_range container_t = std::vector<branch_t>>
class branch_stack
{
private:
    //!\brief The type of the internal container iterator.
    using stack_iterator_type = std::ranges::iterator_t<container_t>;
    //!\brief The inner memory pool to store the branches.
    container_t _stack{};
    //!\brief Iterator to the current statck position.
    stack_iterator_type _stack_iter{_stack.end()};

public:
    /*!\name Associated types
     * \{
     */
    using container_type = container_t; //!< The used container type.
    using value_type = typename container_type::value_type; //!< The value type.
    using size_type = typename container_type::size_type; //!< The size type.
    using reference = typename container_type::reference; //!< The reference type.
    using const_reference = typename container_type::const_reference; //!< The const reference type.
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr branch_stack() = default; //!< Default.
    //!\brief Copy construct from other.
    constexpr branch_stack(branch_stack const & other) : _stack{other._stack}
    {
        _stack_iter = std::ranges::next(_stack.begin(), std::ranges::distance(other._stack.begin(), other._stack_iter));
    }

    //!\brief Moves construt from other.
    constexpr branch_stack(branch_stack && other) noexcept : branch_stack{}
    {
        swap(other);
    }

    //!\brief Copy/move assigns from other.
    constexpr branch_stack & operator=(branch_stack other) noexcept
    {
        swap(other);
        return *this;
    }

    ~branch_stack() = default; //!< Default.
    //!\}

    /*!\brief Element access
     * \{
     */
    //!\brief Returns a reference to the branch at the given position.
    reference branch_at(size_type position) noexcept
    {
        assert(position < size());
        return _stack[position];
    }

    //!\overload
    const_reference branch_at(size_type position) const noexcept
    {
        assert(position < size());
        return _stack[position];
    }

    //!\brief Returns a reference to the base branch.
    reference base_branch() noexcept
    {
        assert(!empty());
        return _stack.front();
    }

    //!\overload
    const_reference base_branch() const noexcept
    {
        assert(!empty());
        return _stack.front();
    }

    //!\brief Returns a reference to the last pushed branch.
    reference top() noexcept
    {
        assert(!empty());
        return *_stack_iter;
    }

    //!\overload
    const_reference top() const noexcept
    {
        assert(!empty());
        return *_stack_iter;
    }
    //!\}

    /*!\name Capacity
     * \{
     */
    bool empty() const noexcept
    {
        return _stack_iter == _stack.end();
    }

    size_type size() const noexcept
    {
        return empty() ? 0 : (std::ranges::distance(_stack.begin(), _stack_iter) + 1);
    }
    //!\}

    /*!\name Modification
     * \{
     */
    //!\brief Prefetches the next item on the stack by reusing already allocated memory.
    reference prefetch() noexcept(noexcept(_stack.resize(1)))
    {
        size_t const old_size = size();
        if (old_size >= _stack.size())
        {
            bool const is_empty = empty();
            _stack.resize(old_size + 1);
            _stack_iter = (is_empty) ? _stack.end() : std::ranges::next(_stack.begin(), old_size - 1);
        }

        return _stack[old_size];
    }

    //!\brief Realises the prefetched element, such that the prefetched element is the next top element of the satck.
    void realise_prefetched() noexcept
    {
        empty() ? _stack_iter = _stack.begin() : ++_stack_iter;
    }

    //!\brief Removes the branch on top of the stack.
    void pop() noexcept
    {
        assert(!empty());
        (_stack_iter == _stack.begin()) ? _stack_iter = _stack.end() : --_stack_iter;
    }

    //!\brief Pushes a new branch on top of the stack.
    void push(branch_t branch) noexcept(noexcept(prefetch()))
    {
        prefetch() = std::move(branch);
        realise_prefetched();
    }

    //!\brief Pushes a new branch on top of the stack using placement construction.
    template <typename ...args_t>
        requires std::constructible_from<branch_t, args_t...>
    void emplace(args_t && ...args) noexcept(noexcept(prefetch()))
    {
        prefetch() = branch_t{std::forward<args_t>(args)...};
        realise_prefetched();
    }

    //!\brief Swaps this with the content of other.
    void swap(branch_stack & other) noexcept
    {
        auto const old_position = std::ranges::distance(other._stack.begin(), other._stack_iter);
        _stack.swap(other._stack);
        _stack_iter = std::ranges::next(_stack.begin(), old_position);
    }
    //!\}
};
}  // namespace libjst::detail
