// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::journaled_sequence_tree_cursor.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <ranges>

#include <libjst/context_position.hpp>

namespace libjst::no_adl
{

//!\brief Specific class implementation in no_adl namespace to avoid ADL of template arguments.
template <typename jst_t>
class journaled_sequence_tree_cursor
{
public:
    class type;
};

//!\brief Implements the actual journaled sequence tree cursor type.
template <typename jst_t>
class journaled_sequence_tree_cursor<jst_t>::type
{
private:

    //!\brief The iterator type over the sequence collection.
    using sequence_collection_iterator = std::ranges::iterator_t<typename jst_t::sequence_collection_type const>;
    //!\brief The iterator type over a particular sequence.
    using sequence_iterator = std::ranges::iterator_t<typename jst_t::sequence_type const>;
    //!\brief The type of the sequence context.
    using sequence_context_type = std::ranges::subrange<sequence_iterator, sequence_iterator>;
    //!\brief The type to store the positions of a particular sequence context.
    using context_positions_type = std::vector<context_position>;

    jst_t const * _jst_ptr; //!< The referenced jst
    size_t _context_size{}; //!< The size of the underlying context.

    sequence_collection_iterator _sequence_collection_it{}; //!< The iterator over the sequence collection.
    sequence_iterator _sequence_window_it_begin{}; //!< The window iterator over a particular sequence.
    sequence_iterator _sequence_window_it_end{}; //!< The window end iterator over a particular sequence.

public:

    using context_position_type = context_position; //!< The type representing a single context position.

    // The iterator.
    class iterator;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr type() = default; //!< Default.
    constexpr type(type const &) = default; //!< Default.
    constexpr type(type &&) = default; //!< Default.
    constexpr type & operator=(type const &) = default; //!< Default.
    constexpr type & operator=(type &&) = default; //!< Default.
    ~type() = default; //!< Default.

    /*!\brief Constructs the journaled sequence tree cursor for a given libjst::journaled_sequence_tree and a context
     *        size.
     *
     * \param[in] jst A pointer to a const journaled sequence tree.
     * \param[in] context_size The context size to use for the cursor.
     *
     * \details
     *
     * The cursor is initialised with the given context size and already represents the first context.
     */
    type(jst_t const * jst, size_t const context_size) noexcept :
        _jst_ptr{jst},
        _context_size{context_size}
    {
        assert(_jst_ptr != nullptr);

        _sequence_collection_it = std::ranges::begin(_jst_ptr->_sequences);
        init_next_sequence_window();
    }
    //!\}

    /*!\name Iterators
     * \{
     */
    //!\brief Returns an iterator for the current context enumerator.
    iterator begin() noexcept
    {
        return iterator{this};
    }

    //!\overload
    iterator begin() const = delete;

    //!\brief Returns a sentinel for the current context enumerator.
    std::default_sentinel_t end()
    {
        return std::default_sentinel;
    }

    //!\overload
    std::default_sentinel_t end() const = delete;
    //!\}

    //!\brief Returns the current sequence context.
    sequence_context_type context() const
    {
        assert(_jst_ptr != nullptr);

        return sequence_context_type{_sequence_window_it_begin, _sequence_window_it_end};
    }

    //!\brief Returns a vector with libjst::context_position's for each sequence that shares this context.
    context_positions_type positions() const
    {
        assert(_jst_ptr != nullptr);

        return {{collection_position(), sequence_position()}};
    }

    //!\brief Checks wether the cursor is at the end of the libjst::journaled_sequence_tree.
    bool at_end() const noexcept
    {
        assert(_jst_ptr != nullptr);

        return _sequence_collection_it == std::ranges::end(_jst_ptr->_sequences);
    }

    //!\brief Advances the cursor to the next valid context.
    void advance() noexcept
    {
        assert(_jst_ptr != nullptr);

        if (_sequence_window_it_end == std::ranges::end(*_sequence_collection_it))
        {
            std::ranges::advance(_sequence_collection_it, 1);
            init_next_sequence_window();
        }
        else
        {
            std::ranges::advance(_sequence_window_it_begin, 1);
            std::ranges::advance(_sequence_window_it_end, 1);
        }
    }
private:

    //!\brief Determins the id if the currently visited sequence.
    size_t collection_position() const noexcept
    {
        return static_cast<size_t>(std::ranges::distance(std::ranges::begin(_jst_ptr->_sequences),
                                                         _sequence_collection_it));
    }

    //!\brief Determins the begin position of the window of the currently visited sequence.
    size_t sequence_position() const noexcept
    {
        return static_cast<size_t>(std::ranges::distance(std::ranges::begin(*_sequence_collection_it),
                                                         _sequence_window_it_begin));
    }

    /*!\brief Initialises the next sequence window.
     *
     * \details
     *
     * This function is called when the cursor moves to the next sequence. It initialises the context window for the
     * next sequence or advances if the sequence is smaller than the window size.
     */
    void init_next_sequence_window() noexcept
    {
        assert(_jst_ptr != nullptr);

        if (_sequence_collection_it == std::ranges::end(_jst_ptr->_sequences))
            return;

        _sequence_window_it_begin = std::ranges::begin(*_sequence_collection_it);
        _sequence_window_it_end = std::ranges::next(_sequence_window_it_begin,
                                                    _context_size,
                                                    std::ranges::end(*_sequence_collection_it));
        if (static_cast<size_t>(std::ranges::distance(_sequence_window_it_begin, _sequence_window_it_end)) <
            _context_size)
        {
            advance();
        }
    }
};

/*!\brief The context iterator.
 * \implements std::input_iterator
 *
 * \details
 *
 * This iterator implements the C++20 input iterator concept. It is not copyable and allows only comparison with
 * the sentinel.
 */
template <typename jst_t>
class journaled_sequence_tree_cursor<jst_t>::type::iterator
{
private:
    //!\brief The pointer to the underlying host.
    type * _host{nullptr};

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The value type.
    using value_type = typename type::sequence_context_type;
    //!\brief The reference type.
    using reference = value_type;
    //!\brief The difference type.
    using difference_type = std::ptrdiff_t;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief The iterator category.
    using iterator_category = std::input_iterator_tag;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    iterator() = default; //!< Default.
    iterator(iterator const &) = delete; //!< Deleted.
    iterator(iterator &&) = default; //!< Default.
    iterator & operator=(iterator const &) = delete; //!< Deleted.
    iterator & operator=(iterator &&) = default; //!< Default.
    ~iterator() = default; //!< Default.

    /*!\brief Constructs an iterator with the host.
     * \param[in] host The pointer to the host.
     */
    explicit iterator(type * host) : _host{host}
    {}
    //!\}

    /*!\name Member access
     * \{
     */
    //!\brief Returns the pointed to context.
    reference operator*() const
    {
        assert(_host != nullptr);

        return _host->context();
    }

    //!\brief Returns the positions of all sequences that share this context.
    typename type::context_positions_type positions() const
    {
        assert(_host != nullptr);

        return _host->positions();
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Increments this by one.
    iterator & operator++() noexcept
    {
        assert(_host != nullptr);

        _host->advance();
        return *this;
    }

    //!\brief Increments this by one.
    void operator++(int) noexcept
    {
        ++(*this);
    }
    //!\}

    /*!\name Comparison
     * \{
     */
    //!\brief Compares this with the sentinel.
    bool operator==(std::default_sentinel_t const &) const noexcept
    {
        assert(_host != nullptr);

        return _host->at_end();
    }
    //!\}
};

} // namespace libjst::no_adl

namespace libjst
{

/*!\brief A cursor over a libjst::journaled_sequence_tree
 *
 * \tparam jst_t The type of the libjst::journaled_sequence_tree.
 *
 * \details
 *
 * This class implements a cursor over the journaled sequence tree.
 * The cursor provides a context interface over the referentially compressed sequences, i.e. contexts between sequences
 * are processed only once. The cursor can be used inside of the search to enumerate all unique sequence contexts.
 */
template <typename jst_t>
using journaled_sequence_tree_cursor = typename libjst::no_adl::journaled_sequence_tree_cursor<jst_t>::type;

}  // namespace libjst
