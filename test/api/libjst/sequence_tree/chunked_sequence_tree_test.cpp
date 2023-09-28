// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>
#include <string>



#include <libjst/sequence_tree/chunked_tree.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/rcms/dna_compressed_multisequence.hpp>
#include <libjst/rcms/rcs_store.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

#include "../mock/rcs_store_mock.hpp"

namespace jst::test::chunked_sequence_tree {

using source_t = std::string;
using variant_t = jst::test::variant<uint32_t, source_t, uint32_t, std::vector<uint32_t>>;

struct expected_root {
    source_t sequence{};
};

struct fixture {
    source_t source{};
    std::vector<variant_t> variants{};
    uint32_t coverage_size{};
    uint32_t chunk_size{};
    uint32_t overlap_size{};
    uint32_t window_size{};
    std::vector<std::vector<source_t>> expected_labels{};

    template <typename stream_t, typename this_t>
        requires std::same_as<std::remove_cvref_t<this_t>, fixture>
    friend stream_t & operator<<(stream_t & stream, this_t &&) {
        stream << "fixture";
        return stream;
    }
};

struct test : public ::testing::TestWithParam<fixture> {
    using coverage_type = libjst::bit_coverage<uint32_t>;
    using coverage_domain_type = libjst::coverage_domain_t<coverage_type>;

    using cms_t = libjst::dna_compressed_multisequence<source_t, coverage_type>;
    using cms_value_t = std::ranges::range_value_t<cms_t>;
    using rcs_store_t = libjst::rcs_store<source_t, cms_t>;
    rcs_store_t _mock;

    void SetUp() override {
        _mock = rcs_store_t{GetParam().source, GetParam().coverage_size};
        coverage_domain_type domain = _mock.variants().coverage_domain();

        std::ranges::for_each(GetParam().variants, [&] (auto var) {
            _mock.add(cms_value_t{libjst::breakpoint{var.position, var.deletion},
                                  var.insertion,
                                  coverage_type{var.coverage, domain}});
        });
    }

    rcs_store_t const & get_mock() const noexcept{
        return _mock;
    }
};

} // namespace jst::test::chunked_sequence_tree

using namespace std::literals;

using fixture = jst::test::chunked_sequence_tree::fixture;
using variant_t = jst::test::chunked_sequence_tree::variant_t;
struct chunked_sequence_tree_test : public jst::test::chunked_sequence_tree::test
{
    using typename jst::test::chunked_sequence_tree::test::rcs_store_t;
    using jst::test::chunked_sequence_tree::test::get_mock;
    using jst::test::chunked_sequence_tree::test::GetParam;


    auto make_forest() const noexcept {
        auto const & rcs_mock = get_mock();
        return rcs_mock | libjst::chunk(GetParam().chunk_size, GetParam().overlap_size)
                        | std::views::transform([&] (auto partial_tree) {
                            return libjst::labelled(std::move(partial_tree)) | libjst::coloured()
                                                                             | libjst::trim(GetParam().window_size)
                                                                             | libjst::merge();
                        });
    }

    template <typename tree_t>
    constexpr void run_test(tree_t tree, size_t chunk_idx) {
        using node_t = libjst::tree_node_t<tree_t>;
        using cargo_t = libjst::tree_label_t<tree_t>;

        auto cvt_to_string = [] (auto seq) -> std::string {
            std::string str;
            for (char c : seq)
                str.push_back(c);
            return str;
        };

        std::vector<std::string> actual_labels{};
        std::stack<node_t> path{};
        path.push(libjst::root(tree));

        std::cout << "Labels: ";
        while (!path.empty()) {
            node_t p = std::move(path.top());
            path.pop();
            cargo_t label = *p;
            if (!std::ranges::empty(label.sequence())) {
                actual_labels.push_back(cvt_to_string(label.sequence()));
                std::cout << actual_labels.back() << " " << std::flush;
            }

            if (auto c_ref = p.next_ref(); c_ref.has_value()) {
                path.push(std::move(*c_ref));
            }
            if (auto c_alt = p.next_alt(); c_alt.has_value()) {
                path.push(std::move(*c_alt));
            }
        }
        std::cout << "\n";

        std::ptrdiff_t expected_count = std::ranges::ssize(GetParam().expected_labels[chunk_idx]);
        std::ptrdiff_t actual_count = std::ranges::ssize(actual_labels);
        EXPECT_EQ(expected_count, actual_count);
        for (std::ptrdiff_t i = 0; i < std::min(expected_count, actual_count); ++i)
            EXPECT_EQ(cvt_to_string(GetParam().expected_labels[chunk_idx][i]), actual_labels[i]) << i;
    }
};

// ----------------------------------------------------------------------------
// Test case definitions
// ----------------------------------------------------------------------------

TEST_P(chunked_sequence_tree_test, traverse) {
    auto forest = make_forest();

    size_t chunk_idx{};
    for (auto && tree : forest) {
        std::string tracepoint = "Failed in chunk: " + std::to_string(chunk_idx);
        SCOPED_TRACE(tracepoint);
        run_test(std::move(tree), chunk_idx);
        // libjst::tree_traverser_base path{tree};
        // auto expected_it = std::ranges::begin(GetParam().expected_labels[chunk_idx]);
        // for (auto it = path.begin(); it != path.end(); ++it) {
        //     auto && label = *it;
        //     if (std::ranges::empty(label.sequence()))
        //         continue;

        //     EXPECT_TRUE(expected_it != std::ranges::end(GetParam().expected_labels[chunk_idx]));
        //     EXPECT_RANGE_EQ(*expected_it, label.sequence());
        //     if (!std::ranges::equal(*expected_it, label.sequence())) {
        //         std::cout << "Chunk idx: " << chunk_idx << " ";
        //         std::cout << "Label pos: " << std::ranges::distance(std::ranges::begin(GetParam().expected_labels[chunk_idx]), expected_it) << "\n";
        //     }
        //     ++expected_it;
        // }
        ++chunk_idx;
    }
    EXPECT_EQ(chunk_idx, std::ranges::size(GetParam().expected_labels));
}

// ----------------------------------------------------------------------------
// Test values
// ----------------------------------------------------------------------------

using namespace std::literals;

INSTANTIATE_TEST_SUITE_P(no_variant_single_chunk, chunked_sequence_tree_test, testing::Values(fixture{
    .source{"AAAAGGGG"s},
    .variants{},
    .coverage_size{4},
    .chunk_size{8},
    .window_size{4},
    .expected_labels{{"AAAAGGGG"s}}
}));

INSTANTIATE_TEST_SUITE_P(no_variant_two_chunks, chunked_sequence_tree_test, testing::Values(fixture{
    .source{"AAAAGGGG"s},
    .variants{},
    .coverage_size{4},
    .chunk_size{4},
    .window_size{2},
    .expected_labels{{"AAAA"s, "GG"s}, {"GGGG"s}}
}));

INSTANTIATE_TEST_SUITE_P(no_variant_three_chunks, chunked_sequence_tree_test, testing::Values(fixture{
    .source{"AAAAGGGG"s},
    .variants{},
    .coverage_size{4},
    .chunk_size{3},
    .window_size{2},
    .expected_labels{{"AAA"s, "AG"s},
                     {"AGG"s, "GG"s},
                     {"GG"s}}
}));

INSTANTIATE_TEST_SUITE_P(two_variants_single_chunk, chunked_sequence_tree_test, testing::Values(fixture{
         //  01234567
    .source{"AAAAGGGG"s},
    .variants{variant_t{.position{1}, .insertion{"C"s}, .deletion{1}, .coverage{0,1}},
              variant_t{.position{5}, .insertion{"T"s}, .deletion{1}, .coverage{0,2}}},
    .coverage_size{4},
    .chunk_size{8},
    .window_size{4},
    .expected_labels{{"A"s,
                       "CAAG"s,
                           "T"s,
                           "G"s,
                       "AAAG"s,
                           "TGG"s,
                           "GGG"s}}
}));

INSTANTIATE_TEST_SUITE_P(two_variants_two_chunks, chunked_sequence_tree_test, testing::Values(fixture{
         //  01234567
    .source{"AAAAGGGG"s},
    .variants{variant_t{.position{1}, .insertion{"C"s}, .deletion{1}, .coverage{0,1}},
              variant_t{.position{5}, .insertion{"T"s}, .deletion{1}, .coverage{0,2}}},
    .coverage_size{4},
    .chunk_size{4},
    .window_size{3},
    .expected_labels{{"A"s,
                        "CAAG"s,
                        "AAA"s,
                           "G"s,
                            "TG"s,
                            "GG"s
                     },
                     {
                        "G"s,
                          "TGG"s,
                          "GGG"s
                     }}
}));

INSTANTIATE_TEST_SUITE_P(two_variants_three_chunks, chunked_sequence_tree_test, testing::Values(fixture{
         //  01234567
    .source{"AAAAGGGG"s},
    .variants{variant_t{.position{1}, .insertion{"C"s}, .deletion{1}, .coverage{0, 1}},
              variant_t{.position{5}, .insertion{"T"s}, .deletion{1}, .coverage{0, 2}}},
    .coverage_size{4},
    .chunk_size{3},
    .window_size{4},
    .expected_labels{{"A"s,
                       "CAAG"s,
                           "T"s,
                           "G"s,
                       "AA"s,
                         "AG"s,
                           "TG"s,
                           "GG"s
                     },
                     {
                      "AG"s,
                        "TGG"s,
                        "G"s,
                         "GG"s
                     },
                     {
                        "GG"s
                     }}
}));

INSTANTIATE_TEST_SUITE_P(two_variants_two_chunks_overlap, chunked_sequence_tree_test, testing::Values(fixture{
         //  01234567
    .source{"AAAAGGGG"s},
    .variants{variant_t{.position{1}, .insertion{"C"s}, .deletion{1}, .coverage{0, 1}},
              variant_t{.position{5}, .insertion{"T"s}, .deletion{1}, .coverage{0, 2}}},
    .coverage_size{4},
    .chunk_size{4},
    .overlap_size{2},
    .window_size{2},
    .expected_labels{{"A"s,
                       "CAA"s,
                       "AAAG"s,
                           "TGG"s,
                           "G"s,
                            "GG"s
                     },
                     {
                        "G"s,
                          "TGG"s,
                          "GGG"s
                     }}
}));
