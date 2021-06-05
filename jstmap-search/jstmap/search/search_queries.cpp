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

    using search_state_t = typename decltype(libjst::shift_or_algorithm{std::declval<pattern_t &>()})::state_type;
    using search_state_manager_t = libjst::search_state_manager_stack<search_state_t>;
    using algorithm_t = decltype(libjst::shift_or_algorithm{std::declval<pattern_t &>(), search_state_manager_t{}});

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
    pattern_verifier(pattern_t pattern) : _pattern{std::move(pattern)}, _algorithm{_pattern}
    {
        _max_step_size = std::ranges::size(_pattern) + _max_error_count;
    }

    template <std::ranges::input_range haystack_t, typename callback_t>
    void operator()(haystack_t && haystack, callback_t && callback)
    {
        if (std::ranges::empty(_pattern))
        {
            auto haystack_it = haystack.begin();
            callback(optimum{}, haystack_it);
        }
        else
        {
            // We call for every step this function
            _algorithm(haystack, [&] (auto const & state, auto const & haystack_it)
            {
                auto & [current_step, best_hit] = _state_manager.state();

                if (_algorithm.verify(state))
                {
                    if (uint16_t error = _algorithm.error_count(); best_hit.error_count > error)
                        best_hit = optimum{.error_count = error, .step_count = current_step};
                }

                if (++current_step == _max_step_size && best_hit.error_count <= _max_error_count)
                    callback(best_hit, haystack_it);
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
std::vector<libjst::context_position> search_queries_(jst_t const & jst,
                                                      seqan::StringSet<raw_sequence_t> const & queries)
{
    assert(!seqan::empty(queries));

    std::vector<libjst::context_position> results{};

    // ----------------------------------------------------------------------------
    // Initialise and run pigeonhole filter
    // ----------------------------------------------------------------------------

    using state_t = typename decltype(libjst::pigeonhole_filter{queries})::state_type;
    libjst::pigeonhole_filter filter{queries, 0.0, libjst::search_state_manager_stack<state_t>{}};

    size_t const fragment_size = filter.qgram_size();

    auto jst_range_agent = jst.range_agent(fragment_size, filter.state_manager()); // already pushing a branch.
    filter(jst_range_agent, [&] (auto const & hit, auto const & haystack_it)
    {
        std::cout << "hit " << hit << " at: " << haystack_it.coordinate() << "\n";
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

        auto jst_range_extender = jst.range_extender(haystack_it.coordinate());

        // We need some mode of verification:
        pattern_verifier suffix_verifier{query_suffix};
        //TODO: What happens if pattern size = 0?
        auto & forward_extender = jst_range_extender.forward_extender(std::ranges::size(query_suffix),
                                                                      suffix_verifier.state_manager());

        suffix_verifier(forward_extender, [&] (auto && best_suffix_hit, auto const & suffix_it)
        {
            // ----------------------------------------------------------------------------
            // Verify query prefix
            // ----------------------------------------------------------------------------

            pattern_verifier prefix_verifier{query_prefix};
            auto & reverse_extender = jst_range_extender.reverse_extender(std::ranges::size(query_prefix),
                                                                          prefix_verifier.state_manager());
            prefix_verifier(reverse_extender, [&] (auto && best_prefix_hit, auto const & prefix_it)
            {
                uint32_t const total_error_count = best_suffix_hit.error_count + best_prefix_hit.error_count;
                seqan3::debug_stream << "\t- match with: " << total_error_count << " errors\n";
                // match_type match{.journal_decorator = prefix_it.context(), // already the full context?
                //                  .total_errors = best_suffix_hit.error_count + best_prefix_hit.error_count,
                //                  .coordinate = suffix_it.coordinate()}; //
                // Duplicates?
                    // - how can we define duplicates
                    // - if they share the same coordinate
            });
        });
    });
    return results;
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
