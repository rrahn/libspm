// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::pigeonhole_filter.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <ranges>

#include <seqan/index.h>

#include <seqan3/alphabet/range/sequence.hpp>
#include <seqan3/range/views/to_rank.hpp>

#include <libjst/search/state_manager_concept.hpp>
#include <libjst/search/state_manager_single.hpp>
#include <libjst/utility/seqan3_to_seqan2_alphabet_adaption.hpp>

namespace libjst
{

template <seqan3::alphabet alphabet_t>
struct pigeonhole_filter_state
{
    uint64_t hash{};
    alphabet_t left_symbol{};
    uint8_t steps{};
};

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
template <typename pattern_collection_t,
          search_state_manager state_manager_t =
            search_state_manager_single<
                pigeonhole_filter_state<std::ranges::range_value_t<std::ranges::range_value_t<pattern_collection_t>>>>>
// //!\cond
//     requires std::ranges::random_access_range<std::ranges::range_value_t<pattern_collection_t>>
// //!\endcond
class pigeonhole_filter
{
private:

    using pattern_t = std::ranges::range_value_t<pattern_collection_t>; //!\brief The pattern type.
    using alphabet_t = std::ranges::range_value_t<pattern_t>; //!\brief The alphabet type.

    using qgram_shape_t = seqan::Shape<alphabet_t, seqan::SimpleShape>;
    using qgram_index_spec_t = seqan::IndexQGram<qgram_shape_t, seqan::OpenAddressing>;
    using qgram_index_t = seqan::Index<pattern_collection_t, qgram_index_spec_t>;

    static constexpr size_t alphabet_size = seqan3::alphabet_size<alphabet_t>; //!\brief The size of the alphabet.

    qgram_index_t _qgram_index{};
    state_manager_t _state_manager{}; //!\brief The state manager to use for the search.
public:

    using state_type = pigeonhole_filter_state<alphabet_t> ; //!< The type of the state.

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr pigeonhole_filter() = default; //!< Default.
    constexpr pigeonhole_filter(pigeonhole_filter const &) = default; //!< Default.
    constexpr pigeonhole_filter(pigeonhole_filter &&) = default; //!< Default.
    constexpr pigeonhole_filter & operator=(pigeonhole_filter const &) = default; //!< Default.
    constexpr pigeonhole_filter & operator=(pigeonhole_filter &&) = default; //!< Default.
    ~pigeonhole_filter() = default; //!< Default.

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
    template <search_state_manager other_state_manager_t = search_state_manager_single<state_type>>
    pigeonhole_filter(pattern_collection_t const & pattern_collection,
                      float const error_rate = 0.0,
                      other_state_manager_t && state_manager = {}) :
        _qgram_index{pattern_collection},
        _state_manager{std::forward<other_state_manager_t>(state_manager)}
    {
        size_t max_length{};
        size_t max_delta = 3;
        size_t min_delta = std::numeric_limits<size_t>::max();
        std::ranges::for_each(pattern_collection, [&] (auto && needle)
        {
            size_t needle_size = std::ranges::size(needle);
            max_length = std::max(max_length, needle_size);
            size_t error_count = static_cast<size_t>(std::floor(error_rate * needle_size));
            size_t delta = needle_size / (error_count + 1);
            if (delta >= 3)
            {
                min_delta = std::min(min_delta, delta);
                max_delta = std::max(max_delta, delta);
            }
        });

        if (min_delta < 3) min_delta = max_delta;

        if (min_delta == std::numeric_limits<size_t>::max()) // disable index
            min_delta = max_length + 1;
        else
            seqan::setStepSize(_qgram_index, min_delta);  // !Important step to change the step size.

        seqan::indexShape(_qgram_index) = qgram_shape_t{std::min<unsigned>(63/3, min_delta)};
        seqan::indexRequire(_qgram_index, seqan::QGramSADir{});


//  Filter initialisation by SeqAn.
//     typedef typename Size<TIndex>::Type TSize;
//     double _newErrorRate = errorRate;
//     TSize seqCount = countSequences(host(pattern));

//     if (pattern._currentErrorRate != _newErrorRate)
//     {
//         // settings have been changed -> initialize bucket parameters

//         pattern._currentErrorRate = _newErrorRate;

//         TSize minDelta = std::numeric_limits<TSize>::max();
//         TSize maxDelta = 3;
//         TSize maxSeqLen = 0;
//         for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
//         {
//             // get pattern length and max. allowed errors
//             TSize length = sequenceLength(seqNo, host(pattern));
//             if (maxSeqLen < length) maxSeqLen = length;

//             // sequence must have sufficient length
//             if (length <= pattern.params.overlap) continue;

//             // cut overlap many characters from the end
//             TSize errors = (TSize) floor(errorRate * length);
//             length -= pattern.params.overlap;
//             TSize delta = length / (errors + 1);


//             // ignore too short q-grams
//             if (delta < 3) continue;
//             if (minDelta > delta) minDelta = delta;
//             if (maxDelta < delta) maxDelta = delta;
//         }
//         pattern.maxSeqLen = maxSeqLen;
//         if (minDelta < 3) minDelta = maxDelta;

//         TIndex &index = host(pattern);
//         pattern.finderPosOffset = 0;
//         pattern.finderPosNextOffset = pattern.maxSeqLen + pattern.finderLength;

//         if (pattern.params.delta != 0)
//         {
//             // use user-defined delta
//             minDelta = pattern.params.delta;
//         }

//         if (minDelta == std::numeric_limits<TSize>::max())
//         {
//             // disable index
//             minDelta = pattern.maxSeqLen + 1;
//         }

//         if (_pigeonholeUpdateShapeLength(pattern.shape, minDelta + pattern.params.overlap) || getStepSize(index) != minDelta)
//         {
//             clear(index);
//             setStepSize(index, minDelta);  !Important step to change the step size.
//          }
//         indexShape(host(pattern)) = pattern.shape;
// //        double start = sysTime();
//         indexRequire(host(pattern), QGramSADir());

        // TODO: construct q_gram index
    }
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
    template <std::ranges::input_range haystack_t,
             typename /*std::invocable<std::ranges::iterator_t<haystack_t> &>*/ callback_t>
    constexpr void operator()(haystack_t && haystack, callback_t && on_hit) noexcept
    {
        // We cannot apply the hash function if we are at the start.
        // But we can already compute the hash value step by step.
        //
        // 01234
        // ####
        // ^####
        // .........
        //     ^
        // update hash: subtract hash value of last element
        // add hash value of current pointer
        // store next character: begin of context.

        qgram_shape_t & shape = seqan::indexShape(_qgram_index);
        size_t const shape_size = seqan::length(shape);

        for (auto haystack_it = std::ranges::begin(haystack); haystack_it != std::ranges::end(haystack); ++haystack_it)
        {
            state_type & state = _state_manager.state();
            if (state.steps >= shape_size) // Compute rolling hash
            {
                // Now we can always run with a rolling hash.
                state.hash = (state.hash - seqan3::to_rank(state.left_symbol) * shape.leftFactor) * alphabet_size;
                state.hash += seqan3::to_rank(*haystack_it);
                state.left_symbol = *(haystack_it.base() - (shape_size - 1));

                process_hash(haystack_it, on_hit);
            }
            else // compute growing hash.
            {
                // If first step ever than set left symbol.
                if (state.steps == 0)
                    state.left_symbol = *haystack_it;

                state.hash = state.hash * alphabet_size + seqan3::to_rank(*haystack_it);

                if (++state.steps == shape_size)
                    process_hash(haystack_it, on_hit);
            }
        }
    }

    //!\brief Returns a reference to the underlying state manager.
    constexpr state_manager_t & state_manager() noexcept
    {
        return _state_manager;
    }

    //!\brief Returns the qgram size.
    constexpr size_t qgram_size() const noexcept
    {
        return seqan::length(seqan::indexShape(_qgram_index));
    }

private:

    template <typename haystack_it, typename callback_t>
    void process_hash(haystack_it const & it, callback_t && on_hit)
    {
        qgram_shape_t & shape = seqan::indexShape(_qgram_index);
        shape.hValue = _state_manager.state().hash;

        for (auto && hit : seqan::getOccurrences(_qgram_index, shape))
            on_hit(hit, it);
    }
};

/*!\name Type deduction guides
 * \relates libjst::pigeonhole_filter
 * \{
 */
//!\brief Deduces the pattern type from the constructor argument.
template <typename collection_t>
pigeonhole_filter(collection_t const &) -> pigeonhole_filter<collection_t>;

//!\brief Deduces the pattern and the state manager type from the constructor arguments.
template <typename collection_t, search_state_manager state_manager_t>
pigeonhole_filter(collection_t const &, float const, state_manager_t &&)
    -> pigeonhole_filter<collection_t, std::decay_t<state_manager_t>>;
//!\}

}  // namespace libjst
