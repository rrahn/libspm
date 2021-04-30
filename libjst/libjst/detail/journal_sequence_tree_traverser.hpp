// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::detail::journal_sequence_tree_traverser.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/zip.hpp>

#include <libjst/context_position.hpp>
#include <libjst/journal_decorator.hpp>
#include <libjst/detail/branch_stack.hpp>
#include <libjst/detail/journal_sequence_tree_traverser_model.hpp>
#include <libjst/utility/logger.hpp>

namespace libjst::detail
{

/*!\brief The basic class to traverse a libjst::journaled_sequence_tree.
 *
 * \tparam derived_t The derived type of this traverser.
 * \tparam jst_t The type of the libjst::journaled_sequence_tree.
 *
 * \details
 *
 * This crtp-base class implements the actual traversal over the libjst::journaled_sequence_tree. It expands each
 * subtree depending on the given context size and the variants contained in the tree. Unsupported branches are not
 * traversed.
 */
template <typename derived_t, typename jst_t>
class journal_sequence_tree_traverser : protected journal_sequence_tree_traverser_model<jst_t>
{
private:
    //!\brief The type of the traverser model.
    using model_t = journal_sequence_tree_traverser_model<jst_t>;

    //!\brief Grant access to the derived class.
    friend derived_t;

    /*!\name Event types
     * \{
     */
    using typename model_t::branch_event_queue_type;
    using typename model_t::branch_event_queue_iterator;
    using typename model_t::branch_event_type;

    using typename model_t::join_event_queue_type;
    using typename model_t::join_event_queue_iterator;
    using typename model_t::join_event_type;
    //!\}

    /*!\name Delta event types
     * \{
     */
    using typename model_t::delta_event_shared_type;
    using segment_type = typename delta_event_shared_type::segment_type; //!< The segment type.
    using typename model_t::coverage_type;
    using typename model_t::size_type;
    using journal_decorator_type = journal_decorator<segment_type>; //!< The journal decorator type.
    //!\}

    //!\brief The type of the sequence context.
    using sequence_context_type = decltype(std::declval<journal_decorator_type const &>() | seqan3::views::slice(0, 1));

    //!\brief The type to store the positions of a particular sequence context.
    using context_positions_type = std::vector<context_position>;

    struct branch;

    //!\brief Represents the status of the branch creation.
    enum branch_creation_status : uint8_t
    {
        success, //!< A new branch could be created.
        no_support, //!< No new branch was created as it was not supported by any of the sequences.
        success_with_deletion //!< A new branch could be created covering a deletion.
    };

    branch_stack<branch> _branch_stack{}; //!< The internal stack of branches.
    join_event_queue_iterator _join_event_it{}; //!< The global join event iterator to update the context positions.
    size_t _context_size{}; //!< The size of the underlying context.

    using context_position_type = context_position; //!< The type representing a single context position.

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr journal_sequence_tree_traverser() = default; //!< Default.
    constexpr journal_sequence_tree_traverser(journal_sequence_tree_traverser const &) = default; //!< Default.
    constexpr journal_sequence_tree_traverser(journal_sequence_tree_traverser &&) = default; //!< Default.
    constexpr journal_sequence_tree_traverser & operator=(journal_sequence_tree_traverser const &)
        = default; //!< Default.
    constexpr journal_sequence_tree_traverser & operator=(journal_sequence_tree_traverser &&) = default; //!< Default.
    ~journal_sequence_tree_traverser() = default; //!< Default.

    /*!\brief Constructs the journaled sequence tree traverser from a given traverser model and a context size.
     *
     * \param[in] model The model to construct the traverser for.
     * \param[in] context_size The context size to use for the cursor.
     *
     * \details
     *
     * Initialises the branch stack by pushing the first branch representing the reference sequence. Then
     * creates the first branches if the current context position is already on a branch event.
     */
    journal_sequence_tree_traverser(model_t model, size_t const context_size) :
        model_t{std::move(model)},
        _context_size{context_size}
    {
        assert(_context_size > 0);
        using insertion_t = typename delta_event_shared_type::insertion_type;

        // ----------------------------------------------------------------------------
        // Prepare the branch and join events for this traverser.
        // ----------------------------------------------------------------------------

        // The following steps initialises the parameters for the context specific traverser.
        size_t const bin_end_position = std::min(this->reference().size(), this->_end_pos + (_context_size - 1));
        _join_event_it = std::ranges::begin(this->join_event_queue());
        branch_event_queue_iterator branch_sentinel = this->branch_event_queue().end();
        if (!this->is_final_bin())
        { // The base branch of this bin will stop at the first non insertion branch at the end position.
            delta_event_shared_type key{bin_end_position, insertion_t{}, coverage_type{}};
            branch_sentinel = this->branch_event_queue().upper_bound(branch_event_type{std::addressof(key)});
        }

        // The next branch event is the first branch if this is the first bin. For all other bins, the context starts
        // fully expanded and the first branch event is the first branch greater or equal than the first context
        // position, i.e. the last position of the context, which is not an insertion if the positions are equal.
        //
        // bin[i - 1]   |         bin[i]                    |           bin[i + 1]
        // _____________[.......:.]_________________________[.......:.]_________________________
        //               --------- <- |context| = 9
        //               ^      ^^
        //      _begin_pos      |last_context_position
        //                      end_pos of bin[i - 1]
        //
        // first_branch_event_it: the first branch event inside of the context which is greater than `_begin_pos`
        //                        and not an insertion, i.e. an insertion at `_begin_pos` is excluded.
        // next_branch_event_it: the first branch event with a branch position greater or equal to
        //                       `last_context_position` and a join position greater than `last_context_position`,
        //                       i.e. insertions at `last_context_position` are excluded.
        //
        // The dotted interval has been fully examined when traversing bin[i - 1] including insertions at
        // `last_context_position`. In bin[i] the traversal starts at _begin_pos but is inside of the base
        // branch as if it had examined all branches in between. Accordingly, all branch events between
        // `first_branch_event_it` and `next_branch_event_it` need to be removed from the base branch coverage.
        branch_event_queue_iterator next_branch_event_it = this->branch_event_queue().begin();
        coverage_type initial_coverage{this->_base_coverage};

        if (!this->is_first_bin())
        {
            // Find the first join event if this is not the first bin.
            delta_event_shared_type join_key{this->_begin_pos, insertion_t{}, coverage_type{}};
            _join_event_it = this->join_event_queue().lower_bound(join_event_type{std::addressof(join_key)});
            size_t last_context_position = std::min(this->_begin_pos + (_context_size - 1), bin_end_position);

            // First event with a branch position >= _begin_pos and a jon position > _begin_pos
            delta_event_shared_type key_first{this->_begin_pos, insertion_t{}, coverage_type{}};
            auto first_branch_event_it =
                this->branch_event_queue().upper_bound(branch_event_type{std::addressof(key_first)});

            // First event with a branch position >= last_context_position and a join position > last_context_position
            delta_event_shared_type key_next{last_context_position, insertion_t{}, coverage_type{}};
            next_branch_event_it = this->branch_event_queue().upper_bound(branch_event_type{std::addressof(key_next)});

            std::ranges::for_each(first_branch_event_it, next_branch_event_it, [&] (branch_event_type const & event)
            {
                initial_coverage.and_not(event.coverage());
            });
        }

        // ----------------------------------------------------------------------------
        // Initialise the base branch covering the reference segment.
        // ----------------------------------------------------------------------------

        _branch_stack.emplace(
            this->_begin_pos,  // current context position.
            bin_end_position, // current branch end position.
            0, // branch offset.
            nullptr, // pointer to the delta event causing the branch.
            next_branch_event_it, // next branch event to consider for this branch.
            branch_sentinel, // The sentinel branch event.
            _join_event_it, // next join event to consider for this branch.
            journal_decorator_type{std::span{this->reference()}}, // the journal decorator of the current branch
            std::move(initial_coverage) // the current branch coverage.
        );
        active_branch().jd_iter = active_branch().journal_decorator.begin();

        // ----------------------------------------------------------------------------
        // Initialise the first branch if any exists at the first position.
        // ----------------------------------------------------------------------------

        while (on_branch_event())
        {
            if (branch_creation_status state = create_branch(); state == branch_creation_status::success)
            { // The branch could be created, i.e. there is a valid branch at position 0 with at least one sequence
              // covering this branch event.
                assert(!is_base_branch());
                assert(active_branch().coverage.any());
                break; // terminate loop since we found a candidate.
            }
            else if (state == branch_creation_status::success_with_deletion)
            { // The branch could be created but is immediately dropped because the deletion is at the beginning.
              // There is no context spanning over the deletion here. The alternative branch has been updated
              // and points to the next branch event.
                assert(!is_base_branch()); // Cannot be in the base branch.
                drop_branch(); // remove branch and continue to see if there is another branch.
            }
        }
    }

    /*!\brief Constructs the journaled sequence tree traverser for a given libjst::journaled_sequence_tree and a context
     *        size.
     *
     * \param[in] jst A pointer to a const journaled sequence tree.
     * \param[in] context_size The context size to use for the cursor.
     * \param[in] begin_pos The begin position of the reference to consider.
     * \param[in] end_pos The end position of the reference to consider.
     *
     * \details
     *
     * Creates first a traverser model from the jst and the reference interval and then constructs the traverser
     * with the given context.
     */
    journal_sequence_tree_traverser(jst_t const * jst,
                                    size_t const context_size,
                                    std::ptrdiff_t begin_pos,
                                    std::ptrdiff_t end_pos) :
        journal_sequence_tree_traverser{model_t{jst, begin_pos, end_pos}, context_size}
    {}
    //!\}

    /*!\name Branch stack operations
     * \{
     */
    //!\brief Returns the active branch.
    branch & active_branch() noexcept
    {
        return _branch_stack.top();
    }

    //!\overload
    branch const & active_branch() const noexcept
    {
        return _branch_stack.top();
    }

    //!\brief Checks if the active branch is the base branch.
    bool is_base_branch() const noexcept
    {
        return _branch_stack.size() == 1u;
    }

    /*!\brief Tests wether the enumerator is at the end.
     *
     * \returns `true`, if no more contexts are available, `false` otherwise.
     */
    bool at_end() const noexcept
    {
        return _branch_stack.empty();
    }

    //!\brief Pushes a new branch on the branch stack.
    void push_branch() noexcept
    {
        _branch_stack.realise_prefetched();
        as_derived()->notify_push();
    }

    //!\brief Removes the current branch from the branch stack.
    void drop_branch() noexcept
    {
        assert(!_branch_stack.empty());

        _branch_stack.pop();
        as_derived()->notify_pop();
    }

    //!\brief Makes a new branch at the current position and switches to this branch.
    branch_creation_status create_branch()
    {
        assert(!at_end());
        assert(active_branch().branch_event_it != std::ranges::end(this->branch_event_queue()));
        assert(on_branch_event());

        // Create a branch from the current one.
        branch & new_branch = _branch_stack.prefetch();
        new_branch = active_branch();

        // Update the delta event and the coverage.
        new_branch.delta_event = active_branch().branch_event_it->event_handle();
        active_branch().branch_event_it = next_branch_event(active_branch());
        update_coverage(new_branch);

        // Terminate early if this branch is not supported by any sequence.
        if (new_branch.coverage.none())
            return branch_creation_status::no_support;

        // Apply the delta event to update the current journal decorator.
        new_branch.branch_event_sentinel = std::ranges::end(this->branch_event_queue());
        record_delta_event(new_branch);
        new_branch.jd_iter = new_branch.journal_decorator.begin();
        if (!new_branch.journal_decorator.empty())
            new_branch.jd_iter += new_branch.context_position;

        new_branch.offset += this->event_offset(new_branch.delta_event);

        // The max end position will be determined and possibly cropped if the natural end of the
        // journal sequence tree is reached before (e.g. the end of the reference sequence).
        size_type const max_end_position = (is_base_branch())
                                         ? new_branch.delta_event->position() + _context_size +
                                           new_branch.delta_event->insertion_size() - 1
                                         : branch_max_end_position();

        new_branch.branch_end_position = std::min(this->max_end_position() + new_branch.offset, max_end_position);
        new_branch.branch_event_it = find_next_relative_branch_event(new_branch);

        // Only return with no support if this branch is a deletion and is at end of the branch because the deletion
        // is at the end of the reference sequence. But try again if there are more events left, because there could
        // be an insertion at the end of the branch such that it is extended with a valid context.
        if (bool is_deletion = new_branch.delta_event->is_deletion(); is_deletion &&
            (new_branch.at_end() && !new_branch.has_more_branch_events()))
        {
            return branch_creation_status::no_support;
        }
        else
        { // Push the new branch onto the stack, which automatically switches the active branch to the new one.
            push_branch();
            return is_deletion ? branch_creation_status::success_with_deletion : branch_creation_status::success;
        }
    }

    /*!\brief Terminates all branches that have been fully visited.
     *
     * \details
     *
     * After advancing the context the end of a branch might have been reached and the consumed branch will be deleted.
     * The next alternative branch might also have been consumed already, e.g. if the previous branch was extended by an
     * insertion close to the end of the reference sequence. Or the alternative branch has no support by any of the
     * sequences, i.e. the coverage is 0 for each sequence.
     */
    void terminate_consumed_branches() noexcept
    {
        assert(!at_end());

        // Skip if branch not at end.
        if (!active_branch().at_end())
            return;

        // Only terminate a branch if it reached the end and either the branch end position is also the branch max
        // end position or there are no more branch events left.
        auto reached_branch_end = [&] () -> bool
        {
            return active_branch().at_end() &&
                   (active_branch().branch_end_position == branch_max_end_position() ||
                   !active_branch().has_more_branch_events());
        };

        // Drop all visited or unsupported branches except the base branch.
        // Check for 0 coverage only in the very last.
        while (!is_base_branch() && (reached_branch_end() || active_branch().coverage.none()))
            drop_branch();
    }
    //!\}

    /*!\name Branch operations
     * \{
     */

    //!\brief Returns the iterator to the next branch event not the end of the branch queue.
    branch_event_queue_iterator next_branch_event(branch const & current_branch) const noexcept
    {
        return std::ranges::next(current_branch.branch_event_it, 1, current_branch.branch_event_sentinel);
    }

    /*!\brief Finds the next relative valid branch event.
     *
     * \returns The branch iterator to the first branch event that comes relatively after the deleted/replaced region.
     *
     * \details
     *
     * To find the next branch event two cases need to be considered. In the first case the current branch iterator
     * points to a substitution or an deletion. The next iterator will then point to the first branch event that comes
     * after the deleted segment in the reference.
     * In the second case, the branch iterator points to an insertion. In this case, all insertions at the same position
     * are first skipped and then the next valid branch is searched. The rational here is that insertions are joined
     * at the same position where they are spawned, because they "consume" no character from the reference sequence.
     * But one sequence cannot cover two insertions at the same position as this would result in an ambiguous sequence
     * context (depending on which insertion is traversed first). Thus, all insertions at the current position are
     * skipped and then the regular search is conducted.
     */
    model_t::branch_event_queue_iterator find_next_relative_branch_event(branch const & new_branch) const noexcept
    {
        auto is_local_insertion = [&] (auto const & event) constexpr
        {
            return event.event_handle()->is_insertion() && event.position() == new_branch.delta_event->position();
        };

        auto it = std::ranges::find_if_not(next_branch_event(new_branch),
                                           std::ranges::end(this->branch_event_queue()),
                                           is_local_insertion);

        auto event_lies_behind_deletion = [&] (auto const & event) constexpr
        {
            return event.position() >= new_branch.delta_event->position() + new_branch.delta_event->deletion_size();
        };

        return std::ranges::find_if(it, std::ranges::end(this->branch_event_queue()), event_lies_behind_deletion);
    }

    //!\brief Tests whether the context position of the active branch lies on a branch event.
    bool on_branch_event() const noexcept
    {
        // Either there are no branches left or the context postion is equal to the postion of the branch event
        // relative to the position of the active branch.
        return active_branch().has_more_branch_events() &&
               active_branch().context_position == active_branch().relative_branch_event_position();
    }

    /*!\brief Returns the maximal end position of a branch.
     *
     * \returns the maximal end position of the current branch.
     *
     * \details
     *
     * Only one branch will be processed at a time. Thus, the origin of the current branch will be at the second
     * position within the branch stack from the end. The maximal branch position can be computed with the stored
     * branch position as well as the context size and the respective insertion size of the original branch.
     */
    size_type branch_max_end_position() const noexcept
    {
        assert(_branch_stack.size() > 1);

        branch const & branch_origin = _branch_stack.branch_at(1);
        return branch_position() + _context_size + branch_origin.delta_event->insertion_size() - 1;
    }

    //!\brief Return the original branch position of the current branch.
    size_type branch_position() const noexcept
    {
        assert(_branch_stack.size() > 1);

        return _branch_stack.branch_at(1).delta_event->position();
    }

    //!\brief Return the original branch event of the current branch.
    delta_event_shared_type const * original_branch_event() noexcept
    {
        assert(_branch_stack.size() > 1);

        return _branch_stack.branch_at(1).delta_event;
    }

    /*!\brief Tests wether the full context is available in the active branch.
     *
     * \returns `true` if the full context can be dereferenced, `false` otherwise.
     *
     * \details
     *
     * The context end position must be greater than the context size plus the set begin position.
     * In case of a partitioned JST the context could be at the begin of a new bin and thus requires the same
     * initialisation as in the beginning of the journal sequence tree.
     */
    bool has_full_context_in_branch() const noexcept
    {
        return context_end_position() >= (_context_size + this->_begin_pos);
    }

    //!\brief Returns the begin position of the context in the active branch.
    size_type context_begin_position() const noexcept
    {
        return active_branch().context_position + 1 - _context_size;
    }

    //!\brief Returns the end position of the context in the active branch.
    size_type context_end_position() const noexcept
    {
        return active_branch().context_position + 1;
    }

    //!\brief Records the represented delta event in the journal decorator of the new branch.
    void record_delta_event(branch & new_branch)
    {
        using insertion_t = delta_event_shared_type::insertion_type;
        using deletion_t = delta_event_shared_type::deletion_type;

        size_type position = new_branch.delta_event->position() + new_branch.offset;
        journal_decorator_type & jd = new_branch.journal_decorator;
        std::visit([&] (auto & delta_kind)
        {
            seqan3::detail::multi_invocable
            {
                [&] (insertion_t const & e) { jd.record_insertion(position, std::span{e.value()}); },
                [&] (deletion_t const & e) { jd.record_deletion(position, position + e.value()); },
                [&] (auto const & e) { jd.record_substitution(position, std::span{e.value()}); }
            } (delta_kind);
        }, new_branch.delta_event->delta_variant());
    }

    /*!\brief Updates the coverages of the new branch and the parent branch.
     *
     * \details
     *
     * If the new branch comes from the base branch it will be set to the original event coverage, otherwise to the
     * conjunction of parant branch coverage and the branch event coverage. Hence, the new branch is only represented
     * by all sequences that share the current branch and the same branch history.
     * The parent branch coverage is updated accordingly, to represent all sequences without the new branch event.
     */
    void update_coverage(branch & new_branch) noexcept
    {
        new_branch.coverage = new_branch.delta_event->coverage();
        if (!is_base_branch())
            new_branch.coverage &= active_branch().coverage;

        active_branch().coverage.and_not(new_branch.delta_event->coverage());
    }

    /*!\brief Updates the relative context position offset for each sequence.
     *
     * \details
     *
     * When the context advances in any branch, the head position of the respective context might lie behind a join
     * event. In this case, the relative position offsets for all sequences are updated with the respective event
     * offset of the pointed-to delta event, if it is shared by the respective sequence id.
     * The first branch whose context head position fulfills the condition issues the update.
     * This update affects all branches globally and is done only once for every join event.
     */
    void update_relative_sequence_offsets() noexcept
    {
        if (_join_event_it == this->join_event_queue().end() || _join_event_it->position() > context_begin_position())
            return;

        join_event_queue_iterator upper_bound_join_it{};
        if (is_base_branch())
        {
            upper_bound_join_it = this->join_event_queue().upper_bound(context_begin_position());
        }
        else // inside of branch
        {
            using insertion_t = typename delta_event_shared_type::insertion_type;
            delta_event_shared_type event{branch_position(), insertion_t{}, coverage_type{}};
            join_event_type join_key{&event};
            if (original_branch_event()->is_insertion()) // Stop at the insertion
                upper_bound_join_it = this->join_event_queue().lower_bound(join_key);
            else // Stop at the event itself or the insertion that ends at the same position if there is one.
                upper_bound_join_it = this->join_event_queue().upper_bound(join_key);
        }
        // Update the sequence offsets.
        std::ranges::for_each(_join_event_it, upper_bound_join_it, [&] (auto const & join_event)
        {
            this->update_offset_for_event(join_event);
        });

        _join_event_it = upper_bound_join_it;
    }

    /*!\brief Keep the base branch coverage up-to-date with the joined branches.
     *
     * \details
     *
     * When the base branch is advanced and the context head position lies on a join event, the coverage of the base
     * branch will be updated with the coverages of all events pointed-to by the join events on this position.
     */
    void update_base_branch_coverage() noexcept
    {
        auto head_on_join_event = [&] ()
        {
            return active_branch().join_event_it != std::ranges::end(this->join_event_queue()) &&
                   context_begin_position() == active_branch().join_event_it->position();
        };

        if (is_base_branch() && head_on_join_event())
        {

            for (; head_on_join_event(); ++active_branch().join_event_it)
                active_branch().coverage |= active_branch().join_event_it->coverage();

            // In the base branch the coverage needs to be refined. When joining delta events the coverage will
            // be updated accordingly, which might overwrite a particular bit in the coverage that was previously
            // unset by another delta event that lies inside of the current context of the base branch.

            // The branch_event_it points to the first branch event that lies behind the tail of the current context.
            // Go back to find first branch event that does not affect the coverage of the current context.
            auto first_it = std::ranges::find_if(std::reverse_iterator{active_branch().branch_event_it},
                                                 std::reverse_iterator{std::ranges::begin(this->branch_event_queue())},
                                                 [&] (auto const & event)
            {
                return (event.event_handle()->is_insertion() && event.position() == context_begin_position()) ||
                       (event.position() < context_begin_position());
            }).base();

            assert(first_it == std::ranges::end(this->branch_event_queue()) ||
                   first_it->position() >= context_begin_position());

            // Now determine the supported coverage within this context.
            std::ranges::for_each(first_it, active_branch().branch_event_it, [&] (auto const & branch_event)
            {
                active_branch().coverage.and_not(branch_event.coverage());
            });
        }
    }

    /*!\brief Advances to the next position.
     *
     * \returns `true` if not at end, `false` otherwise.
     *
     * \details
     *
     * When called the context of the current branch is advanced by one position.
     * After advancing the context the following sequence of operations is performed if applicable:
     * * terminate all consumed branches
     * * update relative context position of each sequence
     * * update base branch coverage with joined branches
     * * create a new branch
     * * terminate base branch if fully consumed
     */
    bool advance()
    {
        assert(!at_end());

        ++active_branch().context_position;
        ++active_branch().jd_iter;
        terminate_consumed_branches(); // terminate all branches
        update_base_branch_coverage();

        // If on a branch event, try to create a new branch.
        while (on_branch_event())
        {
            if (branch_creation_status status = create_branch(); status == branch_creation_status::success)
                break;
            else if (status == branch_creation_status::no_support)
                terminate_consumed_branches();
        }

        // Either on the base branch or there must be at least one supported sequence.
        assert(is_base_branch() || active_branch().coverage.any());

        // Drop the base branch if fully consumed as well.
        if (is_base_branch() && active_branch().at_end() && !active_branch().has_more_branch_events())
            drop_branch();

        return at_end();
    }

    //!\brief Advances the context by one position.
    //!\returns `true` if there is full context available, otherwise `false`.
    bool next_context()
    {
        return advance() || has_full_context_in_branch();
    }

    /*!\brief Returns the current context.
     *
     * \returns A range representing the current context.
     *
     * \details
     *
     * The returned range is a slice over the generated journal decorator of the current branch.
     */
    sequence_context_type current_context() const noexcept
    {
        return active_branch().journal_decorator
             | seqan3::views::slice(context_begin_position(), context_end_position());
    }

    //!\brief Returns the current context iterator.
    std::ranges::iterator_t<journal_decorator_type> current_iterator() const noexcept
    {
        return active_branch().jd_iter;
    }

    //!\brief Returns the current value pointed to by the traverser.
    //!\returns The current pointed-to value.
    std::ranges::range_value_t<journal_decorator_type> current_value() const noexcept
    {
        return *active_branch().jd_iter;
    }

   /*!\brief Compute the branch coverage that is valid for the current branch.
     *
     * \returns The valid branch coverage for the current context.
     *
     * \details
     *
     * The traversal does not keep a valid coverage all the time. When asking for the positions the iterator
     * first refines the current coverage to produce a supported coverage, effectively deferring these operations until
     * they are requested by the caller.
     */
    coverage_type determine_supported_context_coverage() noexcept
    {
        update_relative_sequence_offsets();

        // Just use coverage in these cases.
        if (this->branch_event_queue().empty() || is_base_branch() || context_begin_position() >= branch_position())
            return active_branch().coverage;

        // On non-base branches the coverage will not be updated when joining delta events.
        // Instead, this update is done here only when requested.
        // To do this, the interval of join events which lie between the head of the context and the original
        // branch event is found, and then the supported coverage is determined.
        auto join_begin = std::ranges::find_if(active_branch().join_event_it,
                                               std::ranges::end(this->join_event_queue()),
                                               [&] (auto const & join_event)
        {
            return join_event.position() > context_begin_position();
        });

        // Find the end of the interval.
        auto join_end = std::ranges::find_if(join_begin,
                                             std::ranges::end(this->join_event_queue()),
                                             [&] (auto const & join_event)
        {
            return (join_event.event_handle() == original_branch_event()) ||
                    join_event.position() > branch_position();
        });

        // Update the join event so in case we need to call context position on the same branch again, we can
        // start from a closer event.
        active_branch().join_event_it = join_begin;

        // Determine the supported coverage.
        coverage_type unsupported = std::accumulate(join_begin,
                                                    join_end,
                                                    coverage_type(model_t::_sequence_offsets.size(), false),
                                                    [] (coverage_type coverage, auto const & event)
        {
            return coverage | event.coverage();
        });

        return active_branch().coverage & ~unsupported;
    }
    //!\}

    //!\brief Cast this to its derived type.
    derived_t * as_derived() noexcept
    {
        return static_cast<derived_t *>(this);
    }

    //!\overload
    derived_t const * as_derived() const noexcept
    {
        return static_cast<derived_t const *>(this);
    }
};

//!\brief The type of a traversal branch.
template <typename derived_t, typename jst_t>
struct journal_sequence_tree_traverser<derived_t, jst_t>::branch
{
    size_type context_position{}; //<! The current tail position of the moving window.
    size_type branch_end_position{}; //!<! The end position of the branch.
    std::make_signed_t<size_type> offset{}; //<! The offset generated by the current branch.
    delta_event_shared_type const * delta_event{}; //!< The pointer to the current delta event.
    branch_event_queue_iterator branch_event_it{}; //!< The iterator pointing to the next branch event.
    branch_event_queue_iterator branch_event_sentinel{}; //!< The iterator pointing to branch event sentinel.
    join_event_queue_iterator join_event_it{}; //!< The iterator pointing to the next join event.
    journal_decorator_type journal_decorator{}; //!< The journal decorator representing the current sequence context.
    coverage_type coverage{}; //!< The coverage for this branch.
    std::ranges::iterator_t<journal_decorator_type> jd_iter{}; //!< Iterator points to the current journal decorator.

    //!\brief Checks wether this branch reached its end.
    constexpr bool at_end() const noexcept
    {
        return context_position == branch_end_position;
    }

    //!\brief Returns the position of the delta event.
    constexpr size_type delta_event_position() const noexcept
    {
        return delta_event->position();
    }

    //!\brief Returns the position of the pointed-to branch event.
    constexpr size_type branch_event_position() const noexcept
    {
        return branch_event_it->position();
    }

    //!\brief Returns the position of the pointed-to join event.
    constexpr size_type join_event_position() const noexcept
    {
        return join_event_it->position();
    }

    //!\brief Returns the relative branch event position.
    constexpr size_type relative_branch_event_position() const noexcept
    {
        return branch_event_position() + offset;
    }

    //!\brief Tests wether the branch has more branch events left.
    constexpr bool has_more_branch_events() const noexcept
    {
        return branch_event_it != branch_event_sentinel;
    }
};
}  // namespace libjst::detail
