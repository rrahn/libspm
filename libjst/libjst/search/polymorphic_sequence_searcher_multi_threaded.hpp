// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides polymorphic sequence searcher.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <cmath>

#include <omp.h>

#include <libcontrib/matcher/concept.hpp>
#include <libjst/sequence_tree/chunked_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/traversal/state_capture_traverser.hpp>
#include <libjst/traversal/state_oblivious_traverser.hpp>

namespace libjst
{

    template <typename polymorphic_sequence_t>
    class polymorphic_sequence_searcher_multi_threaded
    {
    private:

        using sequence_type = typename polymorphic_sequence_t::source_type;

        std::reference_wrapper<polymorphic_sequence_t const> _polymorphic_sequence;
        uint32_t _thread_count{1};

    public:
        polymorphic_sequence_searcher_multi_threaded() = delete;

        constexpr explicit polymorphic_sequence_searcher_multi_threaded(
                polymorphic_sequence_t const & polymorphic_sequence,
                uint32_t const thread_count) :
            _polymorphic_sequence{std::cref(polymorphic_sequence)},
            _thread_count{std::max<uint32_t>(1, thread_count)}
        {}

        template <libjst::window_matcher pattern_t, typename callback_t>
            // requires libjst::online_matcher_for<pattern_t, sequence_type const &, callback_t>
        constexpr void operator()(pattern_t && pattern, callback_t && callback) const {

            // std::cout << "thread count = " << _thread_count <<"\n";
            uint32_t chunk_size = std::ceil(std::ranges::size(_polymorphic_sequence.get().source())/_thread_count);
            auto forest = libjst::volatile_tree{_polymorphic_sequence.get()} | chunk(chunk_size);
            // std::cout << "forest size = " << std::ranges::size(forest) <<"\n";

            #pragma omp parallel for num_threads(_thread_count), firstprivate(pattern, callback, _thread_count), shared(forest, _polymorphic_sequence), schedule(auto)
            for (uint32_t i = 0; i < _thread_count; ++i) {
                // std::cout << "thread: " << omp_get_thread_num() << " works on " << omp_get_thread_num() << "\n";
                auto traverser = make_traverser<pattern_t>(_polymorphic_sequence.get());
                traverser(forest[omp_get_thread_num()], pattern, callback);
            }
        }

    private:

        template <typename pattern_t>
        static constexpr auto make_traverser([[maybe_unused]] polymorphic_sequence_t const & polymorphic_sequence) noexcept {
            if constexpr (libjst::restorable_matcher<pattern_t>) {
                if constexpr (libjst::reducable_state<libjst::matcher_state_t<pattern_t>>) {
                    return state_oblivious_traverser{};
                } else {
                    return state_capture_traverser{};
                }
            } else {
                return state_oblivious_traverser{};
            }
        }
    };

}  // namespace libjst
