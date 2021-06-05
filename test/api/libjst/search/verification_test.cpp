// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <libjst/journaled_sequence_tree.hpp>
#include <libjst/search/pigeonhole_filter.hpp>
#include <libjst/search/state_manager_stack.hpp>


#include "../test_utility.hpp"

using namespace std::literals;

using sequence_t = seqan3::dna5_vector;

struct verification_test : public ::testing::Test
{};

// TEST(verification_test, verify_use_case)
// {
//     using seqan3::operator""_dna5;

//     using jst_t = libjst::journaled_sequence_tree<sequence_t>;

//     seqan::StringSet<sequence_t> pattern_collection{};
//     seqan::appendValue(pattern_collection, "acgtaacgtaacgtagacga"_dna5);
//     seqan::appendValue(pattern_collection, "acgtacgactacgtacgact"_dna5);
//     seqan::appendValue(pattern_collection, "acgtacgactagcgactacg"_dna5);

//     std::filesystem::path jst_file{DATADIR"sim_refx5.jst"};
//     jst_t jst = libjst::test::load_jst<jst_t>(jst_file);

//     using state_t = typename decltype(libjst::pigeonhole_filter{pattern_collection})::state_type;
//     // filter set error rate.
//     // set state manager from outside.
//     libjst::pigeonhole_filter filter{pattern_collection, libjst::search_state_manager_stack<state_t>{}};
//     auto jst_range_agent = jst.range_agent(filter.qgram_size(), filter.state_manager()); // already pushing a branch.
//     // sequence_t haystack{"acgtaacgtaacgtagacgaacgtacgactacgtacgactacgtacgactagcgactacg"_dna5};

//     filter(jst_range_agent, [&] (auto const & hit, auto const & haystack_it)
//     {
//         // Now the filter gave us some region:
//         // From hit we can infer the position within the read.
//         // From coordinate we can get the exact context position where the read was found.
//         // Assume normal case: How would we verify
//         auto [pattern_idx, pattern_position] = hit;

//         int64_t reference_position = std::ranges::distance(std::ranges::begin(haystack), haystack_it);
//         int64_t reference_hit_end = reference_position + filter.qgram_size();

//         int64_t pattern_prefix_size = pattern_position;
//         int64_t pattern_suffix_size = pattern_collection[pattern_idx].size() - (pattern_position + filter.qgram_size());
//         // Also include the errors, when verifying with indels.
//         auto left_reference_window = haystack
//                                    | seqan3::views::slice(reference_position - pattern_prefix_size, reference_position);
//         auto right_reference_window = haystack
//                                     | seqan3::views::slice(reference_hit_end, reference_hit_end + pattern_suffix_size);

//         error_count = extend_left(left_reference_window | std::views::reverse,
//                                   pattern_collection[pattern_idx] | seqan3::views::slice(0, pattern_prefix_size)
//                                                                   | std::views::reverse,
//                                   error_count); // Upper bound of error count
//         // We would expect the output of the position including the number of errors.
//         // update the error_count
//         size_t pattern_slice_begin = pattern_position + filter.qgram_size();
//         error_count = extend_right(right_reference_window,
//                                    pattern_collection[pattern_idx] |
//                                    seqan3::views::slice(pattern_slice_begin, pattern_slice_begin +  pattern_suffix_size),
//                                    error_count);
//         if (error_count < max_error)
//             // True positive: -> make match
//         std::cout << "hit: " << hit << " at " << haystack_it.coordinate() << "\n";
//     });
// }

// TEST(verification_test, verify_with_jst_extension)
// {
//     using seqan3::operator""_dna5;

//     using jst_t = libjst::journaled_sequence_tree<sequence_t>;

//     // seqan::StringSet<sequence_t> pattern_collection{};
//     // seqan::appendValue(pattern_collection, "acgtaacgtaacgtagacga"_dna5);
//     // seqan::appendValue(pattern_collection, "acgtacgactacgtacgact"_dna5);
//     // seqan::appendValue(pattern_collection, "acgtacgactagcgactacg"_dna5);

//     std::filesystem::path jst_file{DATADIR"sim_refx5.jst"};
//     jst_t jst = libjst::test::load_jst<jst_t>(jst_file);

//     // using state_t = typename decltype(libjst::pigeonhole_filter{pattern_collection})::state_type;
//     // filter set error rate.
//     // set state manager from outside.
//     // libjst::pigeonhole_filter filter{pattern_collection, libjst::search_state_manager_stack<state_t>{}};
//     // auto jst_range_agent = jst.range_agent(filter.qgram_size(), filter.state_manager()); // already pushing a branch.
//     // sequence_t haystack{"acgtaacgtaacgtagacgaacgtacgactacgtacgactacgtacgactagcgactacg"_dna5};

//     // To test this we could simply create a range traverser:
//     auto jst_range_agent = jst.range_agent(4ull);

//     for (auto it = jst_range_agent.begin(); it != jst_range_agent.end(); ++it)
//     {
//         [[maybe_unused]] auto extender = jst.range_extender(it.coordinate());
//         auto & forward_extender = extender.forward_extender(5);
//         std::cout << "Extending coordinate: " << it.coordinate()  << "\n";
//         for (auto ot = forward_extender.begin(); ot != forward_extender.end(); ++ot)
//         {
//             seqan3::debug_stream << *ot;
//         }
//         std::cout << "\n";
//     }




//     // filter(jst_range_agent, [&] (auto const & hit, auto const & haystack_it)
//     // {
//     //     // Now the filter gave us some region:
//     //     // From hit we can infer the position within the read.
//     //     // From coordinate we can get the exact context position where the read was found.
//     //     // Assume normal case: How would we verify
//     //     auto [pattern_idx, pattern_position] = hit;

//     //     // Needs to be able to set its own stack element.
//     //     // Is one level above the traverser.
//     //     //
//     //     auto extender = jst.extender_agent(haystack.coordinate()); // Do we want to generate a new agent every time we have a hit?
//     //     // so extender gets a coordinate
//     //     // needs to have a traverser
//     //     // set the traverser to the seek position
//     //     // now the forward_extender inherits from this object and allows to set some additional parameter,
//     //         // initialised to the first position after the context -> first returned value, plus might already put stuff on the stack
//     //     // the left traverser just does the same -> it does not "recognises the stack" on creation and uses the principle
//     //     // of reverse iterator. -> So points to current iterator but returns the value before,
//     //     // might be not to fast but could be neglectable as well.
//     //     // could be using it's own stack, reference is always present.

//     //     // if found the error region then we can always just use the
//     //     // we know the reference position
//     //     // we can get two coordinates from this reference position!
//     //     // That allows to potentially retrieve the state from the search.

//     //     // extender now points to the region in the tree, where the hit was found.

//     //     // Now we have a multiplicity: For every branch node during right extension we need to do a left extension

//     //     // Forward extender provides a range which allows to simply provide a sub tree view over the right side of the
//     //     // pattern.
//     //     // The original extender tracks the original start position and ensures, that the pattern is not extended beyond
//     //     // the scope of the initial hit.
//     //     shift_or_pattern_searcher verifier_right{pattern_suffix, max_error_count};
//     //     auto forward_extender = extender.forward_extender(suffix_size, verifier_right.state_manager());
//     //     verifier_right(forward_extender, [] (auto && search_result)
//     //     {
//     //         search_result.haystack_it; // Can give me a position and a context if I want.
//     //         search_result.error_count;

//     //         shift_or_pattern_searcher verifier_left(pattern_prefix_reverse, max_error_count - error_count);
//     //         auto reverse_extender = extender.reverse_extender(prefix_size, verifer_left.state_manager());
//     //         verifier_left(reverse_extender, [] (auto && search_result)
//     //         {
//     //             // Found a match!
//     //             // total errors = error_count_right + error_count_left;
//     //             // slice = journal_decorator | slice
//     //             // postion: (how can we locate this position again?)
//     //             // we could store a negative value inside the coordinate
//     //                 // this means: jump to the event and move x steps left.
//     //                 // or if we know the event we are in: we can go from right there
//     //                 // but then we never now how many steps we need to do.
//     //                 // maybe an interval: with begin/end reference
//     //                 // then we can always go back into this position.
//     //         });
//     //     });

//     //     for (auto forward_it = forward_extender.begin(); forward_it != forward_extender.end(); ++forward_it)
//     //     {
//     //         auto current_state = state_manager.state();
//     //         // Accordingly, I can run the normal pattern matching algorithm, which returns only true when it finds a
//     //         // hit.
//     //         if (current_state.steps == 0 && current_state.error_count < total_error)
//     //         {
//     //             size_t max_error_count = total_error - current_state.error_count;
//     //             auto reverse_extender = extender.reverse_extender(context_size); // initialised at the beginning of the initial coordinate
//     //             for (auto reverse_it = reverse_extender.begin(); extender_it != reverse_extender.end(); ++reverse_it)
//     //             {
//     //                 auto current_state = reverse_state_manager.state();
//     //                 if (current_state.steps == 0 && current_state.error_count < total_error)
//     //                 {
//     //                     on_match(reference_context, number_of_errors); // So we can do a
//     //                 }
//     //                 else
//     //                 {
//     //                     verification_kernel(*reverse_it, *state.pattern_it); // Who stores the verification column?
//     //                     ++state.pattern_it; // reverse_pattern iterator
//     //                 }
//     //             }
//     //         }
//     //         else
//     //         {
//     //             verification_kernel(*state.pattern_it, *forward_it);
//     //             ++state.pattern_it;
//     //         }
//     //     }

//     //     int64_t reference_position = std::ranges::distance(std::ranges::begin(haystack), haystack_it);
//     //     int64_t reference_hit_end = reference_position + filter.qgram_size();

//     //     int64_t pattern_prefix_size = pattern_position;
//     //     int64_t pattern_suffix_size = pattern_collection[pattern_idx].size() - (pattern_position + filter.qgram_size());
//     //     // Also include the errors, when verifying with indels.
//     //     auto left_reference_window = haystack
//     //                                | seqan3::views::slice(reference_position - pattern_prefix_size, reference_position);
//     //     auto right_reference_window = haystack
//     //                                 | seqan3::views::slice(reference_hit_end, reference_hit_end + pattern_suffix_size);

//     //     error_count = extend_left(left_reference_window | std::views::reverse,
//     //                               pattern_collection[pattern_idx] | seqan3::views::slice(0, pattern_prefix_size)
//     //                                                               | std::views::reverse,
//     //                               error_count); // Upper bound of error count
//     //     // We would expect the output of the position including the number of errors.
//     //     // update the error_count
//     //     size_t pattern_slice_begin = pattern_position + filter.qgram_size();
//     //     error_count = extend_right(right_reference_window,
//     //                                pattern_collection[pattern_idx] |
//     //                                seqan3::views::slice(pattern_slice_begin, pattern_slice_begin +  pattern_suffix_size),
//     //                                error_count);
//     //     if (error_count < max_error)
//     //         // True positive: -> make match
//     //     std::cout << "hit: " << hit << " at " << haystack_it.coordinate() << "\n";
//     // });
// }

TEST(verification_test, verify_left)
{
    using seqan3::operator""_dna5;

    using jst_t = libjst::journaled_sequence_tree<sequence_t>;
    using event_t = typename jst_t::event_type;
    using snp_t = typename event_t::snp_type;
    using substitution_t = typename event_t::substitution_type;
    using coverage_t = typename event_t::coverage_type;

    //                01234567890123456789
    sequence_t ref = "acgtacgtacgtacgtacgt"_dna5;
    //                  aa
    //                    c
    jst_t jst{std::move(ref), 4};

    jst.insert(event_t{2ull, substitution_t{"AA"_dna5}, coverage_t{1, 0, 1, 0}});
    jst.insert(event_t{4ull, snp_t{"C"_dna5}, coverage_t{1, 1, 0, 0}});
    jst.insert(event_t{4ull, snp_t{"T"_dna5}, coverage_t{0, 0, 0, 1}});

    // To test this we could simply create a range traverser:
    auto enumerator = jst.context_enumerator(4ull);

    for (auto it = enumerator.begin(); it != enumerator.end(); ++it)
    {
        seqan3::debug_stream << "Extending coordinate: " << it.coordinate() << " context: " << *it <<  "\n";
        auto extender = jst.range_extender(it.coordinate());
        auto & reverse_extender = extender.reverse_extender(5);
        for (auto ot = reverse_extender.begin(); ot != reverse_extender.end(); ++ot)
        {
            seqan3::debug_stream << *ot;
        }
        std::cout << "\n";
    }

}

