// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <cereal/archives/binary.hpp>
#include <numeric>

#include <libjst/search/naive_search.hpp>
#include <jstmap/search/load_queries.hpp>
#include <jstmap/search/search_queries.hpp>

#include <libcontrib/seqan/horspool_pattern.hpp>
#include <libjst/concept.hpp>
#include <libjst/search/concept.hpp>
#include <libjst/search/search_base.hpp>
#include <libjst/journaled_sequence_tree/journaled_sequence_tree_model.hpp>
#include <libjst/journaled_sequence_tree/journaled_sequence_tree_forward.hpp>
#include <libjst/journaled_sequence_tree/serialiser_concept.hpp>
#include <libjst/journaled_sequence_tree/serialiser_direct.hpp>
#include <libjst/journaled_sequence_tree/serialiser_delegate.hpp>
#include <libjst/traversal/searcher_factory.hpp>

#include "benchmark_utility.hpp"

template <typename ...args_t>
static void naive_search_benchmark(benchmark::State & state, args_t && ...args)
{
    auto [haplotypes_file] = std::tuple{args...};

    auto haplotypes = jstmap::load_queries(haplotypes_file);
    sequence_t needle = generate_query(state.range(0));
    size_t const total_bytes = std::ranges::size(needle);

    size_t hit_count{};
    std::vector<std::pair<size_t, int32_t>> hits{};

    for (auto _ : state)
    {
        hits.clear();
        libjst::naive_pattern_searcher searcher{needle};

        for (size_t i = 0; i < 100; ++i)
        {
            for (size_t j = 0; j < haplotypes.size(); ++j)
            {
                auto const & haystack = haplotypes[j];
                searcher(haystack, [&] (auto const & haystack_it)
                {
                    hits.emplace_back(i * j, std::ranges::distance(haystack.begin(), haystack_it));
                });
            }
        }
        hit_count += hits.size();
    }

    benchmark::DoNotOptimize(hit_count);
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(total_bytes);
    state.counters["#hits"] = hit_count;
}

// template <typename ...args_t>
// static void jst_search_benchmark(benchmark::State & state, args_t && ...args)
// {
//     auto [jst_file] = std::tuple{args...};
//     auto jst = jstmap::load_jst(jst_file);

//     jstmap::partitioned_jst_t pjst{std::addressof(jst)};

//     std::vector<sequence_t> query{generate_query(state.range(0))};

//     size_t hit_count{};
//     for (auto _ : state)
//         hit_count += jstmap::search_queries(pjst, query).size();

//     benchmark::DoNotOptimize(hit_count);
//     state.counters["bytes_per_second"] =  seqan3::test::bytes_per_second(jst.total_symbol_count());
//     state.counters["#hits"] = hit_count;
// }

struct collect_hits
{
    size_t &hits;

    template <typename finder_t>
    void set_next(finder_t const &) noexcept
    {
        ++hits;
    }

    void set_value() const noexcept
    {}
};

template <typename ...args_t>
static void jst_search_benchmark2(benchmark::State & state, args_t && ...args)
{
    auto [jst_file] = std::tuple{args...};

    jstmap::raw_sequence_t reference{};
    jstmap::jst_model_t jst_model{reference, 0};
    jstmap::fwd_jst_t jst{jst_model};

    std::ifstream archive_stream{jst_file, std::ios_base::binary | std::ios_base::in};
    {
        cereal::BinaryInputArchive in_archive{archive_stream};
        auto jst_archive = in_archive | libjst::direct_serialiser(reference)
                                      | libjst::delegate_serialiser(jst_model);
        libjst::load(jst, jst_archive);
    }

    // libjst::journaled_sequence_tree fwd{std::move(jst)};

    std::vector<jstmap::raw_sequence_t> query{sample_query(reference, state.range(0))};
    jst::contrib::horspool_pattern pattern{query[0]};
    // auto pattern_state = pattern.search_operation();

    size_t hit_count{};
    for (auto _ : state) {
        libjst::search_base(jst, libjst::jst_searcher(pattern.search_operation()), [&] (auto const &) { ++hit_count; });
    }

    size_t pangenome_bytes = std::ranges::size(reference);
    auto && store = libjst::variant_store(jst);
    pangenome_bytes = std::accumulate(store.begin(), store.end(), pangenome_bytes, [] (size_t bytes, auto const & variant)
    {
        return bytes + std::ranges::size(libjst::insertion(variant));
    });

    state.counters["bytes"] = pangenome_bytes;
    state.counters["bytes_per_second"] =  seqan3::test::bytes_per_second(pangenome_bytes);
    state.counters["#hits"] = hit_count;
}

// Register the function as a benchmark
// BENCHMARK_CAPTURE(naive_search_benchmark,
//                   vcf_indel_test,
//                   DATADIR"sim_ref_10Kb_SNP_INDELs_haplotypes.fasta.gz")->Arg(64)->Arg(100)->Arg(150);

BENCHMARK_CAPTURE(jst_search_benchmark2,
                  vcf_indel_test,
                  "/home/rahn/workspace/jstmap/build/data/1KGP.chr22.vcf_new3.jst")->Arg(16)->Arg(32)->Arg(64)->Arg(128)->Arg(256);
                //   "/home/rahn/workspace/jstmap/build/bench/Release/test1.jst")->Arg(16)->Arg(32)->Arg(64)->Arg(128)->Arg(256);
                //   DATADIR"1KGP.chr22.vcf_new.jst")->Arg(64);
                //   "/home/rahn/workspace/data/jstmap/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.jst")->Arg(64);

BENCHMARK_MAIN();
