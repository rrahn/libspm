// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::myers_pattern_searcher.
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

struct myers_large_state
{
    using word_type = uint64_t;
    // TODO: Use stack allocator?
    std::vector<word_type> VP{};
    std::vector<word_type> VN{};
    word_type scoreMask{};            // the mask with a bit set at the position of the last active cell
    uint32_t lastBlock{};              // the block containing the last active cell
};

struct myers_small_state
{
    using word_type = uint64_t;

    myers_large_state * large_state{};

    word_type VP0{}; // VP[0] (saves one dereferentiation)
    word_type VN0{}; // VN[0]
    uint32_t errors{}; // the current number of errors

    myers_small_state() = default;
    myers_small_state(myers_small_state const & other) :
        VP0{other.VP0},
        VN0{other.VN0},
        errors{other.errors}
    {
        if (other.large_state)
            large_state = new myers_large_state{*other.large_state};
    }

    myers_small_state(myers_small_state && other) noexcept : myers_small_state{}
    {
        swap(other);
    }

    myers_small_state & operator=(myers_small_state other)
    {
        swap(other);
        return *this;
    }

    ~myers_small_state()
    {
        if (large_state)
            delete large_state;
    }

    void swap(myers_small_state & rhs)
    {
        using std::swap;
        swap(large_state, rhs.large_state);
        swap(VP0, rhs.VP0);
        swap(VN0, rhs.VN0);
        swap(errors, rhs.errors);
    }
};

template <seqan3::sequence pattern_t,
          search_state_manager state_manager_t = search_state_manager_single<myers_small_state>,
          bool is_global_alignment = false> // By default we use approximate pattern search
//!\cond
    requires std::ranges::view<pattern_t> && std::ranges::random_access_range<pattern_t>
//!\endcond
class myers_algorithm
{
private:

    using word_type = uint64_t;
    using alphabet_t = std::ranges::range_value_t<pattern_t>; //!\brief The alphabet type.

    // static constexpr word_type WORD_INDEX_HIGH_BIT = seqan::BitsPerValue<word_type>::VALUE - 1;
    static constexpr size_t alphabet_size = seqan3::alphabet_size<alphabet_t>; //!\brief The size of the alphabet.
    static constexpr size_t MACHINE_WORD_SIZE = sizeof(word_type) * 8;

    pattern_t _pattern{}; //!\brief The pattern to search for.
    state_manager_t _state_manager{}; //!\brief The state manager to use for the search.

    std::vector<word_type> bitMasks{}; // encode the needle with bitmasks for each alphabet character
    word_type lastBit{};
    uint32_t _maxErrors{};
    uint32_t blockCount{};
    bool isLongNeedle{false};

public:

    using state_type = myers_small_state; //!< The type of the state.

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr myers_algorithm() = default; //!< Default.
    constexpr myers_algorithm(myers_algorithm const &) = default; //!< Default.
    constexpr myers_algorithm(myers_algorithm &&) = default; //!< Default.
    constexpr myers_algorithm & operator=(myers_algorithm const &) = default; //!< Default.
    constexpr myers_algorithm & operator=(myers_algorithm &&) = default; //!< Default.
    ~myers_algorithm() = default; //!< Default.

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
        requires (!std::same_as<std::decay_t<sequence_t>, myers_algorithm>) && seqan3::sequence<sequence_t>
    //!\endcond
    explicit myers_algorithm(sequence_t && pattern, uint32_t max_errors = 0) :
        _pattern{std::views::all(std::forward<sequence_t>(pattern))},
        _maxErrors(max_errors)
    {
        // initialise state.
        size_t pattern_size = std::ranges::size(_pattern);

        // We assume blockCount of 1 at the moment
        blockCount = (pattern_size + MACHINE_WORD_SIZE - 1) / MACHINE_WORD_SIZE;

        bitMasks.resize((alphabet_size + 1) * blockCount);

        // encoding the letters as bit-vectors
        for (size_t j = 0; j < pattern_size; ++j)
        {
            bitMasks[blockCount * seqan3::to_rank(_pattern[j]) + j / MACHINE_WORD_SIZE] |=
                static_cast<word_type>(1) << (j % MACHINE_WORD_SIZE);
        }

        state_type initial_state{};
        initial_state.VP0 = ~static_cast<word_type>(0);
        initial_state.VN0 = 0;
        initial_state.errors = static_cast<uint32_t>(pattern_size);
        lastBit = static_cast<word_type>(1) << (pattern_size - 1); // This is bullshit if the pattern is empty.

        if (blockCount > 1)
        {
            isLongNeedle = true;
            initial_state.large_state = new myers_large_state{};
            size_t localMaxErrors = std::min<size_t>(_maxErrors, pattern_size - 1);
            initial_state.errors = localMaxErrors + 1;
            initial_state.large_state->scoreMask = (static_cast<word_type>(1) << (localMaxErrors % MACHINE_WORD_SIZE));
            initial_state.large_state->lastBlock = localMaxErrors / MACHINE_WORD_SIZE;
            assert(initial_state.large_state->lastBlock < blockCount);
            initial_state.large_state->VP.resize(blockCount, initial_state.VP0);
            initial_state.large_state->VN.resize(blockCount, initial_state.VN0);
            lastBit = static_cast<word_type>(1) << ((pattern_size + MACHINE_WORD_SIZE - 1) % MACHINE_WORD_SIZE);
        }

        _state_manager.state() = std::move(initial_state);
    }

    template <std::ranges::viewable_range sequence_t, search_state_manager other_state_manager_t>
    //!\cond
        requires (!std::same_as<std::decay_t<sequence_t>, myers_algorithm>) && seqan3::sequence<sequence_t>

    //!\endcond
    myers_algorithm(sequence_t && pattern, uint32_t max_error, other_state_manager_t &&) :
        myers_algorithm{std::forward<sequence_t>(pattern), max_error}
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
        // size_t count = 0;
        for (auto it = std::ranges::begin(haystack); it != std::ranges::end(haystack); ++it)
        {
            if (isLongNeedle)
                runLongNeedle(it);
            else
                runShortNeedle(it);

            callback(_state_manager.state(), it);
        }
    }

    constexpr bool verify(state_type const & state) const noexcept
    {
        if (!isLongNeedle)
            return state.errors <= _maxErrors;

        assert(state.large_state != nullptr);
        myers_large_state const & large_state = *state.large_state;
        return ((large_state.lastBlock == (blockCount - 1)) && (large_state.scoreMask == lastBit));
    }

    constexpr uint32_t error_count() const noexcept
    {
        return state_manager().state().errors;
    }

    //!\brief Returns a reference to the underlying state manager.
    constexpr state_manager_t & state_manager() noexcept
    {
        return _state_manager;
    }

    constexpr state_manager_t const & state_manager() const noexcept
    {
        return _state_manager;
    }

private:
    // ----------------------------------------------------------------------------
    // Function runLongNeedle()
    // ----------------------------------------------------------------------------

    template <typename TIterator>
    void runLongNeedle(TIterator const & hystkIt) noexcept
    {
        myers_small_state & state = state_manager().state();
        assert(state.large_state != nullptr);

        myers_large_state & large_state = *state.large_state;

        word_type X{};
        word_type D0{};
        word_type HN{};
        word_type HP{};
        word_type carryD0{};
        word_type carryHP{};
        word_type carryHN{};

        carryD0 = carryHN = 0;
        carryHP = static_cast<word_type>(is_global_alignment); // FIXME: replace Noting with TSpec

        // if the active cell is the last of it's block, one additional block has to be calculated
        // What, why is this necessary?
        // limit = s.lastBlock + (s.scoreMask >> (64 - 1)) // if last bit of score mask is set, then limit is lastBlock + 1 otherwise lastBlock + 0
        size_t limit = large_state.lastBlock + static_cast<size_t>(large_state.scoreMask >> (MACHINE_WORD_SIZE - 1));
        //     1 << (localMaxErrors % 64) // localMaxErrors can be bigger than wordsize so the last non full score mask
        //     large_state.scoreMask = (TWord)1 << (localMaxErrors % pattern.MACHINE_WORD_SIZE);

        // Limit is used as the indicator for Ukkonnen?
        limit -= (limit == blockCount); // really necessary?
        // if (limit == largePattern.blockCount) // if limit == lastBlock we reduce limit by one?
        //     limit--;

        size_t shift = blockCount * seqan3::to_rank(*hystkIt); // Initial position for the consecutive bit mask vector for each character?
        // [a0, a1, a2, a3, c0, c1, c2, c3, g0, g1, g2, g3, t0, t1, t2, t3] // -> blockCount == 4

        // computing the necessary blocks, carries between blocks following one another are stored
        for (word_type currentBlock = 0; currentBlock <= limit; ++currentBlock)
        {
            X = bitMasks[shift + currentBlock] | large_state.VN[currentBlock];

            word_type temp = large_state.VP[currentBlock] + (X & large_state.VP[currentBlock]) + carryD0;
            if (carryD0 != static_cast<word_type>(0))
                carryD0 = temp <= large_state.VP[currentBlock];
            else
                carryD0 = temp < large_state.VP[currentBlock];

            D0 = (temp ^ large_state.VP[currentBlock]) | X;
            HN = large_state.VP[currentBlock] & D0;
            HP = large_state.VN[currentBlock] | ~(large_state.VP[currentBlock] | D0);

            X = (HP << 1) | carryHP;
            carryHP = HP >> (MACHINE_WORD_SIZE - 1);

            large_state.VN[currentBlock] = X & D0;

            temp = (HN << 1) | carryHN;
            carryHN = HN >> (MACHINE_WORD_SIZE - 1);

            large_state.VP[currentBlock] = temp | ~(X | D0);

            // if the current block is the one containing the last active cell
            // the new number of errors is computed
            // Only update error of last block that contains the maximum
            if (currentBlock == large_state.lastBlock)
            {
                if ((HP & large_state.scoreMask) != static_cast<word_type>(0))
                    state.errors++;
                else if ((HN & large_state.scoreMask) != static_cast<word_type>(0))
                    state.errors--;
            }
        }

        // updating the last active cell
        while (state.errors > _maxErrors)
        {
            if ((large_state.VP[large_state.lastBlock] & large_state.scoreMask) != static_cast<word_type>(0))
                state.errors--;
            else if ((large_state.VN[large_state.lastBlock] & large_state.scoreMask) != static_cast<word_type>(0))
                state.errors++;

            large_state.scoreMask >>= 1;
            if (large_state.scoreMask == static_cast<word_type>(0))
            {
                large_state.lastBlock--;
                if constexpr (is_global_alignment)
                {
                    if (large_state.lastBlock == static_cast<uint32_t>(-1))
                        break;
                }

                large_state.scoreMask = static_cast<word_type>(1) << (MACHINE_WORD_SIZE - 1);
            }
        }

        if (!((large_state.scoreMask == lastBit) && (large_state.lastBlock == (blockCount - 1))))
        {
            large_state.scoreMask <<= 1;
            if (!large_state.scoreMask)
            {
                large_state.scoreMask = 1;
                large_state.lastBlock++;
            }

            if ((large_state.VP[large_state.lastBlock] & large_state.scoreMask) != static_cast<word_type>(0))
                state.errors++;
            else if ((large_state.VN[large_state.lastBlock] & large_state.scoreMask) != static_cast<word_type>(0))
                state.errors--;
        }
    }

    // ----------------------------------------------------------------------------
    // Function runShortNeedle()
    // ----------------------------------------------------------------------------

    template <typename TIterator>
    void runShortNeedle(TIterator const & hystkIt) noexcept
    {
        word_type X{};
        word_type D0{};
        word_type HN{};
        word_type HP{};
        // get the current state here?
        state_type & state = _state_manager.state();

        // computing the blocks
        X = bitMasks[seqan3::to_rank(*hystkIt)] | state.VN0;

        D0 = ((state.VP0 + (X & state.VP0)) ^ state.VP0) | X;
        HN = state.VP0 & D0;
        HP = state.VN0 | ~(state.VP0 | D0);
        X = (HP << 1) | static_cast<word_type>(is_global_alignment); // FIXME: replace Nothing by TSpec
        state.VN0 = X & D0;
        state.VP0 = (HN << 1) | ~(X | D0);

        if ((HP & lastBit) != static_cast<word_type>(0))
            state.errors++;
        else if ((HN & lastBit) != static_cast<word_type>(0))
            state.errors--;

        // return std::pair<size_t, bool>(1, s.errors <= _maxErrors);
    }
};

/*!\name Type deduction guides
 * \relates libjst::myers_algorithm
 * \{
 */
//!\brief Deduces the pattern type from the constructor argument.
template <std::ranges::viewable_range sequence_t>
myers_algorithm(sequence_t &&) -> myers_algorithm<std::views::all_t<sequence_t>>;

//!\brief Deduces the pattern and the state manager type from the constructor arguments.
template <std::ranges::viewable_range sequence_t, search_state_manager state_manager_t>
myers_algorithm(sequence_t &&, uint32_t, state_manager_t &&)
    -> myers_algorithm<std::views::all_t<sequence_t>, std::decay_t<state_manager_t>>;
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
          search_state_manager state_manager_t = search_state_manager_single<myers_small_state>>
//!\cond
    requires std::ranges::view<pattern_t> && std::ranges::random_access_range<pattern_t>
//!\endcond
class myers_pattern_searcher
{
private:
    myers_algorithm<pattern_t, state_manager_t> _algorithm{};
public:

    using state_type = myers_small_state; //!< The type of the state.

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr myers_pattern_searcher() = default; //!< Default.
    constexpr myers_pattern_searcher(myers_pattern_searcher const &) = default; //!< Default.
    constexpr myers_pattern_searcher(myers_pattern_searcher &&) = default; //!< Default.
    constexpr myers_pattern_searcher & operator=(myers_pattern_searcher const &) = default; //!< Default.
    constexpr myers_pattern_searcher & operator=(myers_pattern_searcher &&) = default; //!< Default.
    ~myers_pattern_searcher() = default; //!< Default.

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
        requires (!std::same_as<std::decay_t<sequence_t>, myers_pattern_searcher>) && seqan3::sequence<sequence_t>

    //!\endcond
    myers_pattern_searcher(sequence_t && pattern, uint32_t max_errors = 0, other_state_manager_t && state_manager = {}) :
        _algorithm{std::views::all(std::forward<sequence_t>(pattern)),
                   max_errors,
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
 * \relates libjst::myers_pattern_searcher
 * \{
 */
//!\brief Deduces the pattern type from the constructor argument.
template <std::ranges::viewable_range sequence_t>
myers_pattern_searcher(sequence_t &&) -> myers_pattern_searcher<std::views::all_t<sequence_t>>;

//!\brief Deduces the pattern and the state manager type from the constructor arguments.
template <std::ranges::viewable_range sequence_t, search_state_manager state_manager_t>
myers_pattern_searcher(sequence_t &&, uint32_t, state_manager_t &&)
    -> myers_pattern_searcher<std::views::all_t<sequence_t>, std::decay_t<state_manager_t>>;
//!\}

}  // namespace libjst
