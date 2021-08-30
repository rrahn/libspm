// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::detail::journal_sequence_tree_range_extender_agent.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

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
class journal_sequence_tree_range_extender_agent :
    protected journal_sequence_tree_traverser<journal_sequence_tree_range_extender_agent<jst_t>, jst_t>
{
private:
    //!\brief The base traversal type.
    using base_t = journal_sequence_tree_traverser<journal_sequence_tree_range_extender_agent<jst_t>, jst_t>;
    //!\brief The model type.
    using model_t = typename base_t::model_t;

    //!\brief Grant access to the protected and private data members.
    friend base_t;

    // Imported types.
    using typename base_t::size_type;
    using typename base_t::segment_type;
    using typename base_t::branch;
    using typename base_t::delta_event_shared_type;
    using typename base_t::traversal_direction;
    using typename base_t::journal_decorator_type;
    using typename base_t::position_type;

    // The iterator type.
    template <traversal_direction>
    class range_extender;

    journal_sequence_tree_coordinate _coordinate{};
    size_t _original_context_begin_position{};
    delta_event_shared_type const * _original_branch_root{nullptr};

    std::unique_ptr<range_extender<traversal_direction::forward>> _registered_forward_extender{nullptr};
    std::unique_ptr<range_extender<traversal_direction::reverse>> _registered_reverse_extender{nullptr};

public:

    using forward_range_extender_type = range_extender<traversal_direction::forward>;
    using reverse_range_extender_type = range_extender<traversal_direction::reverse>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr journal_sequence_tree_range_extender_agent() = delete; //!< Default.
    constexpr journal_sequence_tree_range_extender_agent(journal_sequence_tree_range_extender_agent const &) = delete;
        //!< Default.
    constexpr journal_sequence_tree_range_extender_agent(journal_sequence_tree_range_extender_agent &&) = default;
        //!< Default.
    constexpr journal_sequence_tree_range_extender_agent & operator=(journal_sequence_tree_range_extender_agent const &)
        = delete; //!< Default.
    constexpr journal_sequence_tree_range_extender_agent & operator=(journal_sequence_tree_range_extender_agent &&)
        = default; //!< Default.
    ~journal_sequence_tree_range_extender_agent() = default; //!< Default.

    /*!\brief Constructs the position agent for a given libjst::journaled_sequence_tree and a context size.
     *
     * \param[in] jst A pointer to a const journaled sequence tree.
     * \param[in] context_size The context size to use for the cursor.
     *
     * \details
     *
     * The range agent is initialised with the given context size and attaches the observer to the stack registry.
     */
    journal_sequence_tree_range_extender_agent(jst_t const * jst,
                                               journal_sequence_tree_coordinate coordinate) noexcept :
        journal_sequence_tree_range_extender_agent{model_t{jst,
                                                           position_type{0u, 0u},
                                                           position_type{0u, std::numeric_limits<size_t>::max()}},
                                                   std::move(coordinate)}
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
    journal_sequence_tree_range_extender_agent(model_t model, journal_sequence_tree_coordinate coordinate) noexcept :
        base_t{std::move(model), 1},
        _coordinate{std::move(coordinate)}
    {
        assert(!_registered_forward_extender);
        assert(!_registered_reverse_extender);
        // Note the base class must be move initialised here and not in the initialiser list.
        // The class stores unique pointers to created forward and reverse extension agents.
        // These agents maintain the list of observer for the tree traversal.
        this->initialise();
        this->seek(_coordinate); // Seek to the correct position.

        assert(this->has_full_context_in_branch());

        _original_branch_root = this->active_branch().subtree_root;
        _original_context_begin_position = this->context_begin_position();
    }
    //!\}

    //!\brief Retrieves the sequence positions from the given coordinate.
    template <search_stack_observer ...observer_t>
    forward_range_extender_type & forward_extender(size_t const extension_size, observer_t & ...observer)
    {
        // TODO: Reset the extender, when called again.
        //   - What about the old one? -> should not be active anymore.
        _registered_forward_extender.reset(new forward_range_extender_type{*this, extension_size, observer...});
        _registered_forward_extender->advance();
        return *_registered_forward_extender;
    }

    template <search_stack_observer ...observer_t>
    reverse_range_extender_type & reverse_extender(size_t const extension_size, observer_t && ...observer)
    {

        // Copy the current branch and put it on top of the branch
        // initialise the new extended branch to the left
        // return a range over the extended range.
        _registered_reverse_extender.reset(new reverse_range_extender_type{*this, extension_size, observer...});
        _registered_reverse_extender->advance();
        return *_registered_reverse_extender;
    }

private:
    //!\brief NOOP function which does nothing.
    void notify_push(traversal_direction dir) noexcept
    {
        switch (dir)
        {
            case traversal_direction::forward:
            {
                if (_registered_forward_extender) _registered_forward_extender->notify_push();
                break;
            }
            case traversal_direction::reverse:
            {
                if (_registered_reverse_extender) _registered_reverse_extender->notify_push();
                break;
            }
        }
    }

    //!\brief NOOP function which does nothing.
    void notify_pop(traversal_direction dir) noexcept
    {
        switch (dir)
        {
            case traversal_direction::forward:
            {
                if (_registered_forward_extender) _registered_forward_extender->notify_pop();
                break;
            }
            case traversal_direction::reverse:
            {
                if (_registered_reverse_extender) _registered_reverse_extender->notify_pop();
                break;
            }
        }
    }
};

template <typename jst_t>
journal_sequence_tree_range_extender_agent(libjst::detail::journal_sequence_tree_traverser_model<jst_t>,
                                           journal_sequence_tree_coordinate const)
                                           -> journal_sequence_tree_range_extender_agent<jst_t>;

template <typename jst_t>
template <typename journal_sequence_tree_range_extender_agent<jst_t>::traversal_direction direction>
class journal_sequence_tree_range_extender_agent<jst_t>::range_extender :
    protected search_stack_notification_registry
{
private:
    using observer_registry_t = search_stack_notification_registry;

    template <typename>
    friend class journal_sequence_tree_range_extender_agent;

    delta_event_shared_type _nil_root{};
    journal_sequence_tree_range_extender_agent & _host;
    size_t _initial_stack_size{};
public:
    range_extender() = delete;
    range_extender(range_extender const &) = delete;
    range_extender(range_extender &&) = delete;
    range_extender & operator=(range_extender const &) = delete;
    range_extender & operator=(range_extender &&) = delete;
    ~range_extender() = default;

    template <typename ...observer_t>
    range_extender(journal_sequence_tree_range_extender_agent & host,
                   size_t const extension_size,
                   observer_t & ...observer) :
        observer_registry_t{observer...},
        _host{host}
    {
        using coverage_t = typename delta_event_shared_type::coverage_type;
        using deletion_t = typename delta_event_shared_type::deletion_type;
        using position_t = typename delta_event_shared_type::position_type;

        // First, duplicate the current branch.
        duplicate_active_branch();
        branch & top_branch = _host.active_branch();
        position_t nil_root_position{.idx = 0, .offset = _host._original_context_begin_position};

        // Initialise the active branch.
        if constexpr (direction == traversal_direction::forward)
        {
            nil_root_position.offset += _host._coordinate.context_size - 1;
            _host._context_size = extension_size + 2;
            top_branch.branch_end_position = std::min<size_t>(_host.max_end_position() + top_branch.offset,
                                                              nil_root_position.offset + extension_size + 1);
            _host._subtree_steps = 0;
            top_branch.jd_iter = top_branch.journal_decorator.begin() + nil_root_position.offset;
        }
        else // reverse direction
        {
            _host._reverse_context_size = extension_size + 1;
            // TODO:
                // - init join event iterator
                // - reset context position
                // - reset branch_end_position
                // - reset context size
                // - advance in reverse direction

            // _host._context_size = extension_size + 2; // maybe we can also work with negative context size?
            // What is the context position?
            // Since the extension can be longer than the found pattern for which we like to extend (Edit distance).
            // Thus, we have to get the correct position.
            // Coordinate stores reference position and node id.
            // So we go to the reference position into the node and then step into the subtree
            // We always have the branch position from the coordiante.
            top_branch.context_position = _host._original_context_begin_position + 1;
            // size_t next_join_position = std::min<size_t>(_host._original_context_begin_position, _host._coordinate.position);

            // Need second context size for reverse traversal.
            if (_host._original_branch_root != nullptr && (_host._original_context_begin_position > _host._coordinate.position))
                top_branch.join_event_it = _host.find_next_relative_branch_event_(top_branch, _host._original_branch_root); //join_event_queue().upper_bound(next_join_position);
            else
                top_branch.join_event_it = _host.join_event_queue().upper_bound(position_t{.offset = _host._original_context_begin_position});

            top_branch.join_event_sentinel = std::ranges::begin(_host.join_event_queue());

            // Need to set the journal iterator.
            // We set our own nil root.

            top_branch.branch_end_position = std::max(0, static_cast<int32_t>(_host._original_context_begin_position) -
                                                         static_cast<int32_t>(extension_size));
            top_branch.next_branch_position = _host.next_branch_position_(top_branch);
            top_branch.jd_iter = top_branch.journal_decorator.begin() + _host._original_context_begin_position;
        }

        // Create dummy root for the extension at the current context position.
        coverage_t coverage{};
        coverage.resize(_host.sequence_count(), true);
        _nil_root = delta_event_shared_type{nil_root_position, deletion_t{0}, std::move(coverage)};
        top_branch.subtree_root = std::addressof(_nil_root);
        // advance();
    }

    class iterator
    {
    private:
        range_extender * _host{nullptr};
    public:
        using value_type = std::ranges::range_value_t<segment_type>; //!< The value type of the iterator.
        using reference = value_type; //!< The reference type which is not assignable.
        using pointer = void;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::input_iterator_tag;

        iterator() = default;
        iterator(iterator const &) = default;
        iterator(iterator &&) = default;
        iterator & operator=(iterator const &) = default;
        iterator & operator=(iterator &&) = default;
        ~iterator() = default;

        explicit iterator(range_extender * host) : _host{host}
        {}

        //\!\brief Returns the underlying iterator of the journal decorator.
        auto base() const
        {
            return _host->_host.current_iterator();
        }

        reference operator*() const noexcept
        {
            assert(_host != nullptr);
            return _host->current_value();
        }

        auto context() const noexcept
            -> decltype(_host->current_context())
        {
            return _host->current_context();
        }

        iterator & operator++() noexcept
        {
            assert(_host != nullptr);
            _host->advance();
            return *this;
        }

        void operator++(int) noexcept
        {
            ++(*this);
        }

        bool operator==(std::default_sentinel_t const &) const
        {
            assert(_host != nullptr);
            return _host->at_end();
        }
    };

    iterator begin() noexcept
    {
        return iterator{this};
    }

    std::default_sentinel_t end() noexcept
    {
        return std::default_sentinel;
    }

    void notify_push()
    {
        observer_registry_t::notify_push();
    }

    void notify_pop()
    {
        observer_registry_t::notify_pop();
    }

private:

    bool at_end() const noexcept
    {
        if constexpr (direction == traversal_direction::forward)
            return _host._branch_stack.size() == _initial_stack_size || _host.at_end();
        else
            return _host._branch_stack.size() == _initial_stack_size || _host.at_end_();
    }

    auto advance()
    {
        if constexpr (direction == traversal_direction::forward)
            _host.advance();
        else
            _host.advance_reverse();
    }

    auto current_value()
    {
        return _host.current_value();
    }

    auto current_context() const noexcept
    {
        size_t begin_position{};
        size_t end_position{};
        int32_t total_context_size = _host._context_size + _host._reverse_context_size +
                                     _host._coordinate.context_size - 3;

        if constexpr (direction == traversal_direction::forward)
        {
            end_position = _host.active_branch().context_position + 1;
            begin_position = end_position - total_context_size;
        }
        else
        {
            // context_position currently points to one after the dereferenced symbol.
            // Can we have negative values though? -> we should ensure always positive values by updating the internal reference position with one.
            assert(_host.active_branch().context_position > 0);
            begin_position = _host.active_branch().context_position - 1;
            end_position = begin_position + total_context_size;
        }

        return std::tuple{_host.active_branch().journal_decorator, begin_position, end_position};
    }

    //!\brief Duplicate active branch.
    void duplicate_active_branch() noexcept
    {
        _initial_stack_size = _host._branch_stack.size();
        _host._branch_stack.push(_host.active_branch());
    }
};

} // namespace libjst::detail
