// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <iostream>
#include <iterator>
#include <tuple>
#include <ranges>

#include <seqan3/core/debug_stream.hpp>

#include <libjst/search/shift_or_search.hpp>
#include <libjst/search/horspool_search.hpp>
#include <libjst/search/myers_search.hpp>
#include <libjst/search/pigeonhole_filter.hpp>
#include <libjst/search/state_manager_stack.hpp>

#include <jstmap/search/write_results.hpp>
#include <jstmap/search/search_queries.hpp>

namespace jstmap
{

template <std::ranges::random_access_range pattern_t>
class pattern_verifier
{
private:

    using search_state_t = typename decltype(libjst::myers_algorithm{std::declval<pattern_t &>()})::state_type;
    using search_state_manager_t = libjst::search_state_manager_stack<search_state_t>;
    using algorithm_t = libjst::myers_algorithm<std::views::all_t<pattern_t &>, search_state_manager_t, true>;
                                                // >(libjst::myers_algorithm{std::declval<pattern_t &>(), 0, search_state_manager_t{}});

    struct optimum
    {
        uint16_t error_count{};
        uint16_t step_count{};
    };

    struct verification_stack_manager
    {
        using state_t = std::pair<uint16_t, optimum>;
        using stack_container_t = std::vector<state_t>;
        using stack_t = std::stack<state_t, stack_container_t>;

        static constexpr uint16_t _infinity{std::numeric_limits<uint16_t>::max()};

        search_state_manager_t & _search_stack{};
        stack_t _verifier_stack{{state_t{0, optimum{.error_count = _infinity}}}};

        constexpr state_t & state() noexcept
        {
            assert(!_verifier_stack.empty());
            return _verifier_stack.top();
        }

        constexpr state_t const & state() const noexcept
        {
            assert(!_verifier_stack.empty());
            return _verifier_stack.top();
        }

        constexpr void on_push()
        {
            _search_stack.on_push();
            _verifier_stack.push(state());
        }

        constexpr void on_pop()
        {
            assert(!_verifier_stack.empty());

            _verifier_stack.pop();
            _search_stack.on_pop();
        }
    };

    pattern_t _pattern{};
    algorithm_t _algorithm{};
    verification_stack_manager _state_manager{._search_stack = _algorithm.state_manager()};
    uint16_t _max_error_count{}; //!< Maximal allowed error count.
    uint16_t _max_step_size{}; //!< Maximal step size to go.

public:

    pattern_verifier() = default;
    pattern_verifier(pattern_t pattern, uint16_t const max_errors = 0) :
        _pattern{std::move(pattern)},
        _algorithm{_pattern, max_errors},
        _max_error_count{max_errors}
    {
        _max_step_size = std::ranges::size(_pattern) + _max_error_count;
        std::cout << "Created a new verifier with max_errors:  " << _max_error_count << " max_step_size = " << _max_step_size << "\n";
    }

    template <std::ranges::input_range haystack_t, typename callback_t>
    void operator()(haystack_t && haystack, callback_t && callback)
    {
        if (std::ranges::empty(_pattern))
        {
            auto haystack_it = haystack.begin();
            // for (; haystack_it != haystack.end(); ++haystack_it);
            callback(optimum{}, haystack_it);
        }
        else
        {
            // We call for every step this function
            _algorithm(haystack, [&] (auto const & state, auto const & haystack_it)
            {
                auto & [current_step, best_hit] = _state_manager.state();
                ++current_step;
                if (_algorithm.verify(state))
                {
                    std::cout << "So it is calling here\n";
                    if (uint16_t error = _algorithm.error_count(); best_hit.error_count > error)
                    {
                        std::cout << "Old errors: " << best_hit.error_count << "\n";
                        best_hit = optimum{.error_count = error, .step_count = current_step};
                        std::cout << "New errors: " << best_hit.error_count << " step count: " << best_hit.step_count << "\n";
                    }
                }

                if (current_step == _max_step_size && best_hit.error_count <= _max_error_count)
                {
                    std::cout << "Oh oh, where are you?\n";
                    callback(best_hit, haystack_it);
                }
            });
        }
    }

    verification_stack_manager & state_manager() noexcept
    {
        return _state_manager;
    }
};

void process_hits(std::vector<libjst::context_position> & results,
                  std::vector<libjst::context_position> const & cursor_positions)
{
    std::ranges::for_each(cursor_positions, [&] (libjst::context_position const & mapping_position)
    {
        results.push_back(mapping_position);
    });
}

// /*TODO:
//  *  * Construct Pigeonhole filter
//  *  * Apply verification method
//  *  * Live happily ever after.
//  */
std::vector<search_match> search_queries_(jst_t const & jst,
                                          seqan::StringSet<raw_sequence_t> const & queries,
                                          float const error_rate)
{
    assert(!seqan::empty(queries));

    std::vector<search_match> matches{};

    // ----------------------------------------------------------------------------
    // Initialise and run pigeonhole filter
    // ----------------------------------------------------------------------------

    using state_t = typename decltype(libjst::pigeonhole_filter{queries})::state_type;
    libjst::pigeonhole_filter filter{queries, error_rate, libjst::search_state_manager_stack<state_t>{}};

    size_t const fragment_size = filter.qgram_size();

    auto jst_range_agent = jst.range_agent(fragment_size, filter.state_manager()); // already pushing a branch.
    filter(jst_range_agent, [&] (auto const & hit, auto const & haystack_it)
    {
        auto jst_coordinate = haystack_it.coordinate();
        std::cout << "hit " << hit << " at: " << jst_coordinate << "\n";
        // Now we have the number of hits generated by a scan with the pigeonhole filter.
        auto [query_idx, query_position] = hit;

        // ----------------------------------------------------------------------------
        // Prepare query prefix and suffix
        // ----------------------------------------------------------------------------

        size_t const suffix_begin_position = query_position + fragment_size;
        assert(suffix_begin_position <= std::ranges::size(queries[query_idx]));

        auto query_prefix = queries[query_idx] | seqan3::views::take(query_position) | std::views::reverse;
        auto query_suffix = queries[query_idx] | seqan3::views::drop(suffix_begin_position);

        // ----------------------------------------------------------------------------
        // Verify query suffix
        // ----------------------------------------------------------------------------

        auto jst_range_extender = jst.range_extender(jst_coordinate);

        // We need some mode of verification:

        uint16_t const max_error_count = (std::floor(error_rate * std::ranges::size(queries[query_idx])));
        uint32_t const suffix_extension_size = (std::ranges::size(query_suffix) + max_error_count) *
                                                static_cast<uint32_t>(!std::ranges::empty(query_suffix));
        pattern_verifier suffix_verifier{query_suffix, max_error_count};
        // seqan3::debug_stream << "\t- suffix: " << query_suffix << "\n";
        // We would need to extend the range of the suffix.
        auto & forward_extender = jst_range_extender.forward_extender(suffix_extension_size,
                                                                      suffix_verifier.state_manager());

        suffix_verifier(forward_extender, [&] (auto && best_suffix_hit, auto const & suffix_it)
        {
            // ----------------------------------------------------------------------------
            // Verify query prefix
            // ----------------------------------------------------------------------------

            uint16_t const remaining_error_count = max_error_count - best_suffix_hit.error_count;
            uint32_t const prefix_extension_size = (std::ranges::size(query_prefix) + remaining_error_count) *
                                                    static_cast<uint32_t>(!std::ranges::empty(query_prefix));

            pattern_verifier prefix_verifier{query_prefix, remaining_error_count};
            auto & reverse_extender = jst_range_extender.reverse_extender(prefix_extension_size,
                                                                          prefix_verifier.state_manager());
            prefix_verifier(reverse_extender, [&] (auto && best_prefix_hit, auto const & prefix_it)
            {
                uint32_t const total_error_count = best_suffix_hit.error_count + best_prefix_hit.error_count;
                std::cout << "hit " << hit << " at: " << jst_coordinate << "\n";
                seqan3::debug_stream << "\t- match with: " << total_error_count << " errors\n";
                // TODO: This would be the simplest thing to do here.
                // Optimally, we already filter here?
                auto tpl = std::ranges::empty(query_prefix) ? suffix_it.context() : prefix_it.context();
                auto && [seq, begin_pos, end_pos] = tpl;

                seqan3::debug_stream << "begin pos = " << begin_pos << " end pos = " << end_pos << "\n";
                // Why are we extracting the positions anyhow?
                assert(suffix_extension_size >= best_suffix_hit.step_count);
                assert(prefix_extension_size >= best_prefix_hit.step_count);
                end_pos -= (suffix_extension_size - best_suffix_hit.step_count);
                begin_pos += (prefix_extension_size - best_prefix_hit.step_count);

                seqan3::debug_stream << "begin pos = " << begin_pos << " end pos = " << end_pos << "\n";

                matches.emplace_back(std::move(seq), begin_pos, end_pos, jst_coordinate, query_idx, total_error_count);

                // Duplicates?
                    // - how can we define duplicates
                    // - if they share the same coordinate
            });
        });
    });
    return matches;
}

std::vector<libjst::context_position> search_queries(partitioned_jst_t const & partitioned_jst,
                                                     std::vector<raw_sequence_t> const & queries)
{
    assert(!queries.empty());

    std::vector<libjst::context_position> results{};

    std::ranges::for_each(queries, [&] (raw_sequence_t const & query)
    {
        // prepare searcher
        using state_t = typename decltype(libjst::horspool_pattern_searcher{query})::state_type;
        libjst::horspool_pattern_searcher searcher{query, libjst::search_state_manager_stack<state_t>{}};

        for (uint32_t index = 0; index < partitioned_jst.bin_count(); ++index)
        {
            auto jst_range_agent = partitioned_jst.range_agent(libjst::context_size{static_cast<uint32_t>(query.size())},
                                                               libjst::bin_index{index},
                                                               searcher.state_manager()); // already pushing a branch.
            searcher(jst_range_agent, [&] (auto & it)
            {
                process_hits(results, partitioned_jst.sequence_positions_at(it.coordinate()));
            });
        }
    });

    return results;
}

}  // namespace jstmap
