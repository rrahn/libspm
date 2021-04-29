// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::detail::journal_sequence_tree_traverser_model.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <iterator>

#include <seqan3/range/views/slice.hpp>

#include <libjst/journal_decorator.hpp>

namespace libjst::detail
{

/*!\brief Provides the data model behind the traverser.
 *
 * \tparam jst_t The type of the libjst::journaled_sequence_tree.
 *
 * \details
 *
 * Provides the data model behind the traverser. It separates the actual model from the traversal operations.
 * The model can be constructed, maintained and serialised separately. Later the model can be used to initialise
 * different context enumerator agents or range agents.
 */
template <typename jst_t>
class journal_sequence_tree_traverser_model
{
protected:

    /*!\name Event types
     * \{
     */
    //!\brief The type of the branch event queue.
    using branch_event_queue_type = decltype(std::declval<jst_t>()._branch_event_queue);
    //!\brief The type of the branch event queue iterator.
    using branch_event_queue_iterator = std::ranges::iterator_t<branch_event_queue_type const>;
    //!\brief The type of the branch event.
    using branch_event_type = std::iter_value_t<branch_event_queue_iterator>;

    //!\brief The type of the join event queue.
    using join_event_queue_type = decltype(std::declval<jst_t>()._join_event_queue);
    //!\brief The type of the join event queue iterator.
    using join_event_queue_iterator = std::ranges::iterator_t<join_event_queue_type const>;
    //!\brief The type of the join event.
    using join_event_type = std::iter_value_t<join_event_queue_iterator>;
    //!\}

    /*!\name Delta event types
     * \{
     */
    using delta_event_shared_type = typename jst_t::delta_event_shared_type; //!< The shared delta event type.
    using coverage_type = typename delta_event_shared_type::coverage_type; //!< The coverage type.
    using size_type = typename delta_event_shared_type::size_type; //!< The size type.
    //!\}

    jst_t const * _jst_host; //!< The referenced jst
    std::vector<int32_t> _sequence_offsets{}; //!< The context position offsets.
    coverage_type _base_coverage{}; //!< The initial base coverage for this model.

    size_t _begin_pos;
    size_t _end_pos;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr journal_sequence_tree_traverser_model() = default; //!< Default.
    constexpr journal_sequence_tree_traverser_model(journal_sequence_tree_traverser_model const &)
        = default; //!< Default.
    constexpr journal_sequence_tree_traverser_model(journal_sequence_tree_traverser_model &&) = default; //!< Default.
    constexpr journal_sequence_tree_traverser_model & operator=(journal_sequence_tree_traverser_model const &)
        = default; //!< Default.
    constexpr journal_sequence_tree_traverser_model & operator=(journal_sequence_tree_traverser_model &&)
        = default; //!< Default.
    ~journal_sequence_tree_traverser_model() = default; //!< Default.

    /*!\brief Constructs the model from the journal sequence tree and an interval over the reference to consider.
     *
     * \param[in] jst A pointer to a const journal sequence tree.
     * \param[in] begin_pos The begin position of the reference sequence.
     * \param[in] end_pos The end position of the reference sequence.
     *
     * \details
     *
     * The model is initialised by setting the sequence offsets to the start of the given begin position.
     * Furthermore, the base coverage will be updated, such that the traversal starting at any given slice will
     * always return only valid contexts positions.
     */
    journal_sequence_tree_traverser_model(jst_t const * jst,
                                          std::ptrdiff_t begin_pos,
                                          std::ptrdiff_t end_pos) noexcept :
        _jst_host{jst}
    {
        assert(_jst_host != nullptr);
        assert(begin_pos < end_pos);

        _begin_pos = std::max<size_t>(0, begin_pos);
        _end_pos = std::min<size_t>(end_pos, _jst_host->reference().size());
        // Initialise the context position offsets and the first join event for updating the context positions.
        _sequence_offsets.resize(_jst_host->size(), 0);

        // Scan over all branch events before the begin position and
        // compute the offset from all events ending before the begin and
        // the coverage from all events ending after the begin.
        _base_coverage.resize(_sequence_offsets.size(), true);
        auto first_candidate_it = branch_event_queue().lower_bound(_begin_pos);
        std::ranges::for_each(std::ranges::begin(branch_event_queue()), first_candidate_it,
        [&] (auto const & branch_event)
        {
            if(branch_event.position() + branch_event.event_handle()->deletion_size() <= _begin_pos)
                update_offset_for_event(branch_event);
            else
                _base_coverage &= branch_event.coverage();
        });
    }
    //!\}

protected:
    //!\brief Returns the event queue containing the branch events from the underlying host.
    branch_event_queue_type const & branch_event_queue() const noexcept
    {
        assert(_jst_host != nullptr);

        return _jst_host->_branch_event_queue;
    }

    //!\brief Returns the event queue containing the join events from the underlying host.
    join_event_queue_type const & join_event_queue() const noexcept
    {
        assert(_jst_host != nullptr);

        return _jst_host->_join_event_queue;
    }

    //!\brief Returns the event queue from the underlying host.
    auto reference() const noexcept
    {
        assert(_jst_host != nullptr);

        return _jst_host->reference() | seqan3::views::slice(_begin_pos, _end_pos);
    }

    /*!\brief Returns the relative delta event offset.
     *
     * returns An integer value representing the relative offset generated by this delta event.
     *
     * \details
     *
     * Depending on the kind of the delta event the returned value can be negative, 0, or positive.
     * A negative value represents the number of deleted characters with respect to the reference sequence, a zero means
     * that with respect to the reference sequence no character was deleted or added and a positive value accounts for
     * the number of characters inserted relative to the reference sequence.
     */
    auto event_offset(delta_event_shared_type const * delta_event) const noexcept
    {
        using signed_size_t = std::make_signed_t<size_type>;

        return static_cast<signed_size_t>(delta_event->insertion_size()) -
               static_cast<signed_size_t>(delta_event->deletion_size());
    }

    /*!\brief Updates the sequence offsets for the given event.
     *
     * \tparam event_t The type of the event (branch event or join event).
     *
     * \param[in] event The event to update the relative sequence positions for.
     *
     * \details
     *
     * Updates the relative sequence positions for every sequence that covers this event. If the event is a
     * substitution it will not invoke the update as the relative offset is not affected.
     */
    template <typename event_t>
    void update_offset_for_event(event_t const & event) noexcept
    {
        if (event.event_handle()->is_substitution())
            return;

        //TODO: Vectorise, mask_add?
        auto const offset = event_offset(event.event_handle());
        for (unsigned idx = 0; idx < _sequence_offsets.size(); ++idx)
            _sequence_offsets[idx] += (event.coverage()[idx] ? offset : 0);
    }
};

}  // namespace libjst::detail
