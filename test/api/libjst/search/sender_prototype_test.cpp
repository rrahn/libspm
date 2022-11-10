// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

// #include <string>

// #include <seqan3/core/debug_stream.hpp>
// #include <seqan3/test/expect_range_eq.hpp>

// #include <libcontrib/seqan/alphabet.hpp>
// #include <libcontrib/seqan/horspool_pattern.hpp>

// #include <libjst/concept.hpp>
// #include <libjst/variant/variant_snp.hpp>
// #include <libjst/variant/variant_store_composite.hpp>
// #include <libjst/variant/variant_store_covered.hpp>
// #include <libjst/container/jst_forward.hpp>
// #include <libjst/container/jst_base.hpp>
// #include <libjst/search/search_base.hpp>
// #include <libjst/traversal/concept_execution.hpp>
// #include <libjst/traversal/searcher_factory.hpp>
// #include <libjst/utility/bit_vector.hpp>

// // #include "../journal_sequence_tree_traversal_test_template.hpp"

// // struct sender_prototype_test : public libjst::test::traversal_fixture_base
// // {};

// struct cout_publisher
// {
//     template <typename this_t, typename result_t>
//         requires std::same_as<std::remove_cvref_t<this_t>, cout_publisher>
//     friend void tag_invoke(std::tag_t<libjst::set_next>, this_t &&, result_t && result) noexcept
//     {
//         auto const & [node, finder] = result;
//         std::cout << "hit at: " << seqan::position(finder) << "\n";
//     }

//     template <typename this_t>
//         requires std::same_as<std::remove_cvref_t<this_t>, cout_publisher>
//     friend void tag_invoke(std::tag_t<libjst::set_value>, this_t &&) noexcept
//     {
//         std::cout << "search finished!\n";
//     }

//     template <typename this_t>
//         requires std::same_as<std::remove_cvref_t<this_t>, cout_publisher>
//     friend void tag_invoke(std::tag_t<libjst::set_done>, this_t &&) noexcept
//     {
//         std::cout << "search aborted!\n";
//     }

//     template <typename this_t, typename error_t>
//         requires std::same_as<std::remove_cvref_t<this_t>, cout_publisher>
//     friend void tag_invoke(std::tag_t<libjst::set_error>, this_t &&, error_t const &) noexcept
//     {
//         std::cout << "search errored!\n";
//     }
// };

// TEST(sender_prototype_test, sender_api)
// {
//     using jst::contrib::operator""_dna4;
//     // we need to define a pattern.
//     using alphabet_t = jst::contrib::dna4;
//     using sequence_t = std::vector<alphabet_t>;
//                         //01234567890123456789012345
//     sequence_t reference{"ACGTGACTAGCATCTAGCATCACGAT"_dna4};

//     jst::contrib::horspool_pattern pattern{"ATCACGAT"_dna4};
//     auto pattern_state = libjst::search_operation_old(pattern);
//     EXPECT_EQ(libjst::window_size(pattern_state), 8u);

//     // create the jst model!
//     using snp_store_t = libjst::variant_store_composite<std::vector<libjst::snp_variant<alphabet_t>>>;
//     using covered_store_t = libjst::variant_store_covered<snp_store_t, libjst::bit_vector<>>;
//     using jst_model_t = libjst::jst_base<sequence_t, covered_store_t>;
//     using fwd_jst_t = libjst::jst_forward<jst_model_t>;

//     jst_model_t jst_model{reference, 4}; // empty basic jst model
//     fwd_jst_t fwd_jst{jst_model}; // now the elements should be sorted! // also if this is empty!

//     libjst::search_base(fwd_jst, libjst::jst_searcher(pattern_state), [] (auto && hit) {
//         auto const & [node, finder] = hit;
//         std::cout << "hit at: " << seqan::position(finder) << "\n";
//     });
//     // auto sender = libjst::jst_searcher(pattern_state);
//     // auto op = libjst::connect(sender, cout_publisher{});
//     // libjst::start(op);

//     // step 1: construct lazy searcher pipeline:
//     // ext::pattern p{"needle"s};
//     // // even independent of each other as long as they can communicate via sane interface
//     // auto searcher = p | jst_search_adapter | filter_pigeonhole | extend_right | extend_left | ...;

//     // // step 2: create text structure:
//     // auto jst = load_fwd_jst();
//     // // jst_search adapter handles this.

//     // // jst.search(searcher)? // how do we tell if this is possible!
//     // // that is defined via the interfaces of the objects.
//     // // so we need some adaption of how we can run the same algorithms on different inputs.

//     // // step 3: search
//     // // just delightful
//     // libjst::search(jst, searcher, [] (auto res) { ... client code ...} )
//     // {
//     //     return lazy_result_range{jst, searcher}; // starts executing when calling begin() -> then uses internal buffer to keep results.
//     //     // or alternatively
//     //     asyn_wait(searcher(jst)); // runs the pipeline in blocking manner
//     //     // auto op = sender.connect(async_wait)
//     // }

// }
