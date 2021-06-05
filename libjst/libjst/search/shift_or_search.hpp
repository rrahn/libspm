// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::shift_or_pattern_searcher.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <ranges>

#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/to_rank.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <libjst/search/state_manager_concept.hpp>
#include <libjst/search/state_manager_single.hpp>

namespace libjst
{

template <seqan3::sequence pattern_t,
          search_state_manager state_manager_t = search_state_manager_single<std::vector<uint64_t>>>
//!\cond
    requires std::ranges::view<pattern_t> && std::ranges::random_access_range<pattern_t>
//!\endcond
class shift_or_algorithm
{
private:

    using alphabet_t = std::ranges::range_value_t<pattern_t>; //!\brief The alphabet type.

    static constexpr size_t alphabet_size = seqan3::alphabet_size<alphabet_t>; //!\brief The size of the alphabet.

    std::vector<uint64_t> _mask_table{}; //!\brief The bit mask table.
    pattern_t _pattern{}; //!\brief The pattern to search for.
    state_manager_t _state_manager{}; //!\brief The state manager to use for the search.
    size_t block_count{};
    uint64_t _hit_mask{};
public:

    using state_type = std::vector<uint64_t>; //!< The type of the state.

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr shift_or_algorithm() = default; //!< Default.
    constexpr shift_or_algorithm(shift_or_algorithm const &) = default; //!< Default.
    constexpr shift_or_algorithm(shift_or_algorithm &&) = default; //!< Default.
    constexpr shift_or_algorithm & operator=(shift_or_algorithm const &) = default; //!< Default.
    constexpr shift_or_algorithm & operator=(shift_or_algorithm &&) = default; //!< Default.
    ~shift_or_algorithm() = default; //!< Default.

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
        requires (!std::same_as<std::decay_t<sequence_t>, shift_or_algorithm>) && seqan3::sequence<sequence_t>
    //!\endcond
    explicit shift_or_algorithm(sequence_t && pattern) :
        _pattern{std::views::all(std::forward<sequence_t>(pattern))}
    {
        // initialise state.
        size_t pattern_size = std::ranges::distance(_pattern);
        _hit_mask = 1 << ((pattern_size - 1) & 63);
        block_count = (pattern_size + 63) / 64;
        _mask_table.resize(block_count * alphabet_size);
        uint64_t mask{1};
        int32_t current_block_idx{-1};
        for (size_t pattern_idx = 0; pattern_idx < pattern_size; ++pattern_idx)
        {
            if (pattern_idx & 63) // no at end of block
            {
                mask <<= 1;
            }
            else // at end of block, so initialise new block.
            {
                mask = 1;
                ++current_block_idx;
            }
            auto rank = seqan3::to_rank(_pattern[pattern_idx]);
            size_t block_position = rank * block_count + current_block_idx;
            _mask_table[block_position] |= mask;
        }

        for (auto & mask : _mask_table)
            mask = ~mask;

        _state_manager.state().resize(block_count, ~0ull);
    }

    template <std::ranges::viewable_range sequence_t, search_state_manager other_state_manager_t>
    //!\cond
        requires (!std::same_as<std::decay_t<sequence_t>, shift_or_algorithm>) && seqan3::sequence<sequence_t>

    //!\endcond
    shift_or_algorithm(sequence_t && pattern, other_state_manager_t &&) :
        shift_or_algorithm{std::forward<sequence_t>(pattern)}
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
        constexpr uint64_t carry_mask = (1ull << 63);
        for (auto it = std::ranges::begin(haystack); it != std::ranges::end(haystack); ++it)
        {
            state_type & search_state = _state_manager.state();

            bool carry_bit{};
            for (size_t block_idx = 0; block_idx < block_count; ++block_idx)
            {
                uint64_t & block = search_state[block_idx];
                bool next_carry_bit = block & carry_mask;
                block = ((block << 1) | carry_bit) |
                        _mask_table[block_count * seqan3::to_rank(*it) + block_idx];
                carry_bit = next_carry_bit;
            }

            callback(search_state.back(), it);
        }
    }

    constexpr bool verify(uint64_t const & state) const noexcept
    {
        return !(state & _hit_mask);
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
 * \relates libjst::shift_or_algorithm
 * \{
 */
//!\brief Deduces the pattern type from the constructor argument.
template <std::ranges::viewable_range sequence_t>
shift_or_algorithm(sequence_t &&) -> shift_or_algorithm<std::views::all_t<sequence_t>>;

//!\brief Deduces the pattern and the state manager type from the constructor arguments.
template <std::ranges::viewable_range sequence_t, search_state_manager state_manager_t>
shift_or_algorithm(sequence_t &&, state_manager_t &&)
    -> shift_or_algorithm<std::views::all_t<sequence_t>, std::decay_t<state_manager_t>>;
//!\}


/*!\brief A shift or pattern search implementation.
 *
 * \tparam pattern_t The type of the pattern to search; must model std::ranges::view and seqan3::sequence.
 * \tparam state_manager_t The type of the state manager; must model libjst::search_state_manager.
 *
 * \details
 *
 * Implements a shift or search algorithm by scanning the text from left to right. It generates
 * a bit mask table to compare the search state with. The state is required to be a single integral value but will
 * be internally resolved to a vector type.
 */
template <seqan3::sequence pattern_t,
          search_state_manager state_manager_t = search_state_manager_single<std::vector<uint64_t>>>
//!\cond
    requires std::ranges::view<pattern_t> && std::ranges::random_access_range<pattern_t>
//!\endcond
class shift_or_pattern_searcher
{
private:
    shift_or_algorithm<pattern_t, state_manager_t> _algorithm{};
public:

    using state_type = std::vector<uint64_t>; //!< The type of the state.

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr shift_or_pattern_searcher() = default; //!< Default.
    constexpr shift_or_pattern_searcher(shift_or_pattern_searcher const &) = default; //!< Default.
    constexpr shift_or_pattern_searcher(shift_or_pattern_searcher &&) = default; //!< Default.
    constexpr shift_or_pattern_searcher & operator=(shift_or_pattern_searcher const &) = default; //!< Default.
    constexpr shift_or_pattern_searcher & operator=(shift_or_pattern_searcher &&) = default; //!< Default.
    ~shift_or_pattern_searcher() = default; //!< Default.

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
        requires (!std::same_as<std::decay_t<sequence_t>, shift_or_pattern_searcher>) && seqan3::sequence<sequence_t>

    //!\endcond
    shift_or_pattern_searcher(sequence_t && pattern, other_state_manager_t && state_manager = {}) :
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
 * \relates libjst::shift_or_pattern_searcher
 * \{
 */
//!\brief Deduces the pattern type from the constructor argument.
template <std::ranges::viewable_range sequence_t>
shift_or_pattern_searcher(sequence_t &&) -> shift_or_pattern_searcher<std::views::all_t<sequence_t>>;

//!\brief Deduces the pattern and the state manager type from the constructor arguments.
template <std::ranges::viewable_range sequence_t, search_state_manager state_manager_t>
shift_or_pattern_searcher(sequence_t &&, state_manager_t &&)
    -> shift_or_pattern_searcher<std::views::all_t<sequence_t>, std::decay_t<state_manager_t>>;
//!\}

}  // namespace libjst
