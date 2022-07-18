// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::horspool_pattern_searcher.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <ranges>

#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/to_rank.hpp>

#include <libjst/search/state_manager_concept.hpp>
#include <libjst/search/state_manager_single.hpp>

namespace libjst
{

/*!\brief A horspool pattern search implementation.
 *
 * \tparam pattern_t The type of the pattern to search; must model std::ranges::view and seqan3::sequence.
 * \tparam state_manager_t The type of the state manager; must model libjst::search_state_manager.
 *
 * \details
 *
 * Implements a horspool search algorithm by scanning the pattern and comparing it from right to left. It generates
 * a occurrence table to jump to a position without checking all locations for a pattern occurrence.
 * The state is required to be a single integral value.
 */
template <seqan3::sequence pattern_t,
          search_state_manager state_manager_t = search_state_manager_single<size_t>>
//!\cond
    requires std::ranges::view<pattern_t> && std::ranges::random_access_range<pattern_t>
//!\endcond
class horspool_algorithm
{
private:

    using alphabet_t = std::ranges::range_value_t<pattern_t>; //!\brief The alphabet type.

    static constexpr size_t alphabet_size = seqan3::alphabet_size<alphabet_t>; //!\brief The size of the alphabet.

    std::array<size_t, alphabet_size> occurrence_table{}; //!\brief The occurrence table.
    pattern_t _pattern{}; //!\brief The pattern to search for.
    state_manager_t _state_manager{}; //!\brief The state manager to use for the search.
public:

    using state_type = size_t; //!< The type of the state.

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr horspool_algorithm() = default; //!< Default.
    constexpr horspool_algorithm(horspool_algorithm const &) = default; //!< Default.
    constexpr horspool_algorithm(horspool_algorithm &&) = default; //!< Default.
    constexpr horspool_algorithm & operator=(horspool_algorithm const &) = default; //!< Default.
    constexpr horspool_algorithm & operator=(horspool_algorithm &&) = default; //!< Default.
    ~horspool_algorithm() = default; //!< Default.

    /*!\brief Constructs a new naive pattern searcher.
     *
     * \tparam sequence_t The type of the pattern; must model std::ranges::viewable_range and seqan3::sequence.
     * \tparam other_state_manager_t The type of the state manager to use for the search.
     *
     * \param[in] pattern The pattern to search.
     * \param[in] state_manager The state manager to use. Defaults to libjst::search_state_manager_single
     *
     * \details
     *
     * Initialises the occurrence table by generating the safe jumps for each symbol in the pattern.
     */
    template <std::ranges::viewable_range sequence_t>
    //!\cond
        requires (!std::same_as<std::decay_t<sequence_t>, horspool_algorithm>) && seqan3::sequence<sequence_t>
    //!\endcond
    explicit horspool_algorithm(sequence_t && pattern) :
        _pattern{std::views::all(std::forward<sequence_t>(pattern))}
    {
        // initialise state.
        size_t max_rank = std::ranges::distance(_pattern) - 1;
        occurrence_table.fill(max_rank);
        size_t symbol_idx{};
        for (auto rank : _pattern | seqan3::views::to_rank | std::views::take(max_rank))
            occurrence_table[rank] = max_rank - ++symbol_idx;

        _state_manager.state() = max_rank;
    }

    template <std::ranges::viewable_range sequence_t, search_state_manager other_state_manager_t>
    //!\cond
        requires (!std::same_as<std::decay_t<sequence_t>, horspool_algorithm>) && seqan3::sequence<sequence_t>
    //!\endcond
    horspool_algorithm(sequence_t && pattern, other_state_manager_t &&) :
        horspool_algorithm{std::forward<sequence_t>(pattern)}
    {}
    //!\}

    /*!\brief Invokes the pattern search on the given input range.
     *
     * \tparam haystack_t The type of the haystack to search; must model std::ranges::input_range.
     * \tparam callback_t The type of the callback to invoke on a hit; must model std::invocable with the iterator
     *                    of the haystack.
     *
     * \param[in] haystack The haystack to search the pattern in.
     * \param[in] on_hit The callback to invoke when a hit is found.
     *
     * \details
     *
     * Iterates over the haystack and checks at every position if the pattern can be found withot any errors.
     * If the pattern was found the callback will be invoked with the iterator pointing to the current position.
     */
    template <std::ranges::input_range haystack_t, typename callback_t>
    constexpr void operator()(haystack_t && haystack, callback_t && callback) noexcept
    {
        size_t const pattern_size = std::ranges::distance(_pattern);
        for (auto it = std::ranges::begin(haystack); it != std::ranges::end(haystack); ++it)
        {
            size_t res = _state_manager.state();
            if (res > 0)  // Advance the haystack iterator until we can compare again.
            {
                --_state_manager.state();
            }
            else
            {
                assert(res == 0);

                auto base_iterator = it.base(); // Returns the original iterator to the JD!

                for (res = 1; res < pattern_size; ++res, --base_iterator) {
                    if (_pattern[pattern_size - res] != *base_iterator)
                        break;
                }
                // Store the current jump for this position.
                _state_manager.state() = occurrence_table[seqan3::to_rank(*it)] + 1;
            }

            callback(res == pattern_size, it);
        }
    }

    constexpr bool verify(bool const & state) const noexcept
    {
        return state;
    }

    constexpr uint32_t error_count() const noexcept
    {
        return 0;
    }

    //!\brief Returns a reference to the underlying state manager.
    constexpr state_manager_t & state_manager() noexcept
    {
        return _state_manager;
    }
};

/*!\name Type deduction guides
 * \relates libjst::horspool_algorithm
 * \{
 */
//!\brief Deduces the pattern type from the constructor argument.
template <std::ranges::viewable_range sequence_t>
horspool_algorithm(sequence_t &&) -> horspool_algorithm<std::views::all_t<sequence_t>>;

//!\brief Deduces the pattern and the state manager type from the constructor arguments.
template <std::ranges::viewable_range sequence_t, search_state_manager state_manager_t>
horspool_algorithm(sequence_t &&, state_manager_t &&)
    -> horspool_algorithm<std::views::all_t<sequence_t>, std::decay_t<state_manager_t>>;
//!\}

template <seqan3::sequence pattern_t,
          search_state_manager state_manager_t = search_state_manager_single<size_t>>
//!\cond
    requires std::ranges::view<pattern_t> && std::ranges::random_access_range<pattern_t>
//!\endcond
class horspool_pattern_searcher
{
private:
    horspool_algorithm<pattern_t, state_manager_t> _algorithm{};
public:

    using state_type = size_t; //!< The type of the state.

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr horspool_pattern_searcher() = default; //!< Default.
    constexpr horspool_pattern_searcher(horspool_pattern_searcher const &) = default; //!< Default.
    constexpr horspool_pattern_searcher(horspool_pattern_searcher &&) = default; //!< Default.
    constexpr horspool_pattern_searcher & operator=(horspool_pattern_searcher const &) = default; //!< Default.
    constexpr horspool_pattern_searcher & operator=(horspool_pattern_searcher &&) = default; //!< Default.
    ~horspool_pattern_searcher() = default; //!< Default.

    /*!\brief Constructs a new naive pattern searcher.
     *
     * \tparam sequence_t The type of the pattern; must model std::ranges::viewable_range and seqan3::sequence.
     * \tparam other_state_manager_t The type of the state manager to use for the search.
     *
     * \param[in] pattern The pattern to search.
     * \param[in] state_manager The state manager to use. Defaults to libjst::search_state_manager_single
     *
     * \details
     *
     * Initialises the occurrence table by generating the safe jumps for each symbol in the pattern.
     */
    template <std::ranges::viewable_range sequence_t,
              search_state_manager other_state_manager_t = search_state_manager_single<state_type>>
    //!\cond
        requires (!std::same_as<std::decay_t<sequence_t>, horspool_pattern_searcher>) && seqan3::sequence<sequence_t>

    //!\endcond
    horspool_pattern_searcher(sequence_t && pattern, other_state_manager_t && state_manager = {}) :
        _algorithm{std::views::all(std::forward<sequence_t>(pattern)),
                   std::forward<other_state_manager_t>(state_manager)}
    {}
    //!\}

    /*!\brief Invokes the pattern search on the given input range.
     *
     * \tparam haystack_t The type of the haystack to search; must model std::ranges::input_range.
     * \tparam callback_t The type of the callback to invoke on a hit; must model std::invocable with the iterator
     *                    of the haystack.
     *
     * \param[in] haystack The haystack to search the pattern in.
     * \param[in] on_hit The callback to invoke when a hit is found.
     *
     * \details
     *
     * Iterates over the haystack and checks at every position if the pattern can be found withot any errors.
     * If the pattern was found the callback will be invoked with the iterator pointing to the current position.
     */
    template <std::ranges::input_range haystack_t, std::invocable<std::ranges::iterator_t<haystack_t> &> callback_t>
    constexpr void operator()(haystack_t && haystack, callback_t && on_hit) noexcept
    {
        _algorithm(haystack, [&] (auto const & state, auto & it)
        {
            if (_algorithm.verify(state))
                on_hit(it);
        });
    }

    //!\brief Returns a reference to the underlying state manager.
    constexpr state_manager_t & state_manager() noexcept
    {
        return _algorithm.state_manager();
    }
};

/*!\name Type deduction guides
 * \relates libjst::horspool_pattern_searcher
 * \{
 */
//!\brief Deduces the pattern type from the constructor argument.
template <std::ranges::viewable_range sequence_t>
horspool_pattern_searcher(sequence_t &&) -> horspool_pattern_searcher<std::views::all_t<sequence_t>>;

//!\brief Deduces the pattern and the state manager type from the constructor arguments.
template <std::ranges::viewable_range sequence_t, search_state_manager state_manager_t>
horspool_pattern_searcher(sequence_t &&, state_manager_t &&)
    -> horspool_pattern_searcher<std::views::all_t<sequence_t>, std::decay_t<state_manager_t>>;
//!\}

}  // namespace libjst
