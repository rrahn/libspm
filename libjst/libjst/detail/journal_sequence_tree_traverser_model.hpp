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

#include <cereal/types/vector.hpp>

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

    //!\brief The jst position type.
    using position_type = typename jst_t::position_type;

    jst_t const * _jst_host{}; //!< The referenced jst
    position_type _begin_pos{}; //!< The begin position of this traverser model.
    position_type _end_pos{}; //!< The end position of this traverser model.

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
                                          position_type begin_pos,
                                          position_type end_pos) noexcept :
        _jst_host{jst},
        _begin_pos{std::move(begin_pos)},
        _end_pos{std::move(end_pos)}
    {
        assert(_jst_host != nullptr);
        assert(_begin_pos < _end_pos);
        assert(_begin_pos.idx == _end_pos.idx);

        _end_pos.offset = std::min<size_t>(_end_pos.offset, max_end_position());
    }
    //!\}

protected:

    //!\brief Returns the begin position of this model.
    constexpr size_t begin_position() const noexcept
    {
        return _begin_pos.offset;
    }

    //!\brief Returns the end position of this model.
    constexpr size_t end_position() const noexcept
    {
        return _end_pos.offset;
    }

    //!\brief Return the contig index.
    constexpr size_t contig_index() const noexcept
    {
        return _begin_pos.idx;
    }

    //!\brief Returns `true` if the end position is equal to the reference size, otherwise `false`.
    bool is_final_bin() const noexcept
    {
        return end_position() == max_end_position();
    }

    //!\brief Returns `true` if the begin position is `0`, otherwise `false`.
    bool is_first_bin() const noexcept
    {
        return begin_position() == 0u;
    }

    //!\brief Returns the maximal end position of the underlying journal sequence tree.
    size_t max_end_position() const noexcept
    {
        assert(contig_index() < _jst_host->reference().size());
        return _jst_host->reference_at(contig_index()).size();
    }

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

    branch_event_queue_iterator contig_branch_begin() const noexcept
    {
        assert(_jst_host != nullptr);

        size_t contig_idx = contig_index();

        if (contig_idx == 0)
            return std::ranges::begin(branch_event_queue());
        else
            return branch_event_queue().upper_bound(position_type{--contig_idx, std::numeric_limits<size_t>::max()});
    }

    branch_event_queue_iterator contig_branch_end() const noexcept
    {
        assert(_jst_host != nullptr);

        size_t contig_idx = contig_index();

        if (contig_idx == _jst_host->reference().size() - 1)
            return std::ranges::end(branch_event_queue());
        else
            return branch_event_queue().lower_bound(position_type{++contig_idx, 0u});
    }

    join_event_queue_iterator contig_join_begin() const noexcept
    {
        assert(_jst_host != nullptr);

        size_t contig_idx = contig_index();

        if (contig_idx == 0)
            return std::ranges::begin(join_event_queue());
        else
            return join_event_queue().upper_bound(position_type{--contig_idx, std::numeric_limits<size_t>::max()});
    }

    join_event_queue_iterator contig_join_end() const noexcept
    {
        assert(_jst_host != nullptr);

        size_t contig_idx = contig_index();

        if (contig_idx == _jst_host->reference().size() - 1)
            return std::ranges::end(join_event_queue());
        else
            return join_event_queue().lower_bound(position_type{++contig_idx, 0u});
    }

    //!\brief Returns the event queue from the underlying host.
    auto const & reference() const noexcept
    {
        assert(_jst_host != nullptr);
        assert(contig_index() < _jst_host->reference().size());

        return _jst_host->reference_at(contig_index());
    }

    //!\brief Returns the number of contained sequences.
    size_type sequence_count() const noexcept
    {
        assert(_jst_host != nullptr);

        return _jst_host->size();
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
     * \param[in,out] offsets The vector with the offsets to update.
     * \param[in] event The event to update the relative sequence positions for.
     *
     * \details
     *
     * Updates the relative sequence positions for every sequence that covers this event. If the event is a
     * substitution it will not invoke the update as the relative offset is not affected.
     */
    template <typename event_t>
    void update_offset_for_event(std::vector<int32_t> & offsets, event_t const & event) const
    {
        if (event.event_handle()->is_substitution())
            return;

        //TODO: Vectorise, mask_add?
        auto const offset = event_offset(event.event_handle());
        for (unsigned idx = 0; idx < offsets.size(); ++idx)
            offsets[idx] += (event.coverage()[idx] ? offset : 0);
    }

public:
    /*!\name Serialisation
     * \{
     */
    /*!\brief Saves this traverser model to the given output archive.
     *
     * \tparam output_archive_t The type of the output_archive; must model seqan3::cereal_output_archive.
     *
     * \param[in, out] archive The archive to serialise this object to.
     */
    template <seqan3::cereal_output_archive output_archive_t>
    void save(output_archive_t & archive) const
    {
        archive(_begin_pos, _end_pos);
    }

    /*!\brief Loads this traverser model from the given input archive.
     *
     * \tparam input_archive_t The type of the input_archive; must model seqan3::cereal_input_archive.
     *
     * \param[in, out] archive The archive to serialise this object from.
     * \param[in] jst A pointer to the associated jst.
     */
    template <seqan3::cereal_input_archive input_archive_t>
    void load(input_archive_t & archive, jst_t const * jst)
    {
        assert(jst != nullptr);
        _jst_host = jst;
        archive(_begin_pos, _end_pos);
    }
    //!\}
};

}  // namespace libjst::detail
