// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::journal_sequence_tree_position_agent.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <numeric>
#include <ranges>

#include <seqan3/range/views/zip.hpp>

#include <libjst/context_position.hpp>
#include <libjst/detail/journal_sequence_tree_traverser.hpp>
#include <libjst/search/stack_notification_registry.hpp>
#include <libjst/utility/logger.hpp>

#include <seqan3/core/debug_stream.hpp>

namespace libjst::detail
{

/*!\brief A range over a libjst::journaled_sequence_tree with an integrated libjst::search_stack_notification_registry.
 * \implements std::ranges::input_range
 *
 * \tparam jst_t The type of the libjst::journaled_sequence_tree.
 *
 * \details
 *
 * This agent provides a range interface to the algorithms. During the traversal a stack is used to keep track of which
 * branch is currently visited inside of the algorithm. To allow external algorithms to keep track of these state
 * changes they can attach a libjst::search_stack_observer during the construction.
 * Subsequently, they will be notified whenever a state changes is applied during the traversal.
 */
template <typename jst_t>
class journal_sequence_tree_position_agent :
    protected journal_sequence_tree_traverser<journal_sequence_tree_position_agent<jst_t>, jst_t>
{
private:
    //!\brief The base traversal type.
    using base_t = journal_sequence_tree_traverser<journal_sequence_tree_position_agent<jst_t>, jst_t>;
    //!\brief The model type.
    using model_t = typename base_t::model_t;

    //!\brief Grant access to the protected and private data members.
    friend base_t;

    // Imported types.
    using typename base_t::coverage_type;
    using typename base_t::size_type;
    using typename base_t::segment_type;
    using typename base_t::traversal_direction;
    using typename base_t::position_type;

    // The iterator type.
    class iterator;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr journal_sequence_tree_position_agent() = default; //!< Default.
    constexpr journal_sequence_tree_position_agent(journal_sequence_tree_position_agent const &) = default;
        //!< Default.
    constexpr journal_sequence_tree_position_agent(journal_sequence_tree_position_agent &&) = default;
        //!< Default.
    constexpr journal_sequence_tree_position_agent & operator=(journal_sequence_tree_position_agent const &)
        = default; //!< Default.
    constexpr journal_sequence_tree_position_agent & operator=(journal_sequence_tree_position_agent &&)
        = default; //!< Default.
    ~journal_sequence_tree_position_agent() = default; //!< Default.

    /*!\brief Constructs the position agent for a given libjst::journaled_sequence_tree and a context size.
     *
     * \param[in] jst A pointer to a const journaled sequence tree.
     * \param[in] context_size The context size to use for the cursor.
     *
     * \details
     *
     * The range agent is initialised with the given context size and attaches the observer to the stack registry.
     */
    journal_sequence_tree_position_agent(jst_t const * jst, size_t const context_size) noexcept :
        journal_sequence_tree_position_agent{model_t{jst,
                                                     position_type{0u, 0u},
                                                     position_type{0u, std::numeric_limits<size_t>::max()}},
                                             context_size}
    {}

    /*!\brief Constructs the range agent from a given traverser model and a context size.
     *
     * \param[in] model The model to construct the traverser for.
     * \param[in] context_size The context size to use for the cursor.
     *
     * \details
     *
     * The position agent is initialised with the given traverser model and context size and attaches the observer
     * to the stack registry.
     */
    journal_sequence_tree_position_agent(model_t model, size_t const context_size) noexcept :
        base_t{std::move(model), context_size}
    {
        this->initialise();
    }
    //!\}

    //!\brief Retrieves the sequence positions from the given coordinate.
    auto retrieve_positions(libjst::journal_sequence_tree_coordinate const & coordinate)
        -> decltype(this->retrieve_positions())
    {
        this->seek(coordinate);
        return base_t::retrieve_positions();
    }

private:
    //!\brief NOOP function which does nothing.
    void notify_push(traversal_direction const &) const noexcept
    {}

    //!\brief NOOP function which does nothing.
    void notify_pop(traversal_direction const &) const noexcept
    {}
};

} // namespace libjst::detail
