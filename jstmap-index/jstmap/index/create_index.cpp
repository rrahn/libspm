// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief The method to create the index.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <seqan3/search/views/kmer_hash.hpp>

#include <libjst/sequence_tree/chunked_tree.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/left_extend_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/prune_unsupported.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

#include <jstmap/index/create_index.hpp>

namespace jstmap
{
// TODO: put functionality into class, so that we can configure it.
seqan3::interleaved_bloom_filter<> create_index(rcs_store_t const & rcs_store, index_options const & options)
{
    // now here we need to change the classes.
    auto forest = rcs_store | libjst::chunk(options.bin_size, options.bin_overlap);

    size_t const bin_count = std::ranges::size(forest);
    size_t ibf_size = 2ull * 1024ull * 1024ull * 1024ull; // 2GiBi
    size_t computed_bin_size = ibf_size / (((bin_count + 63) / 64) * 64);
    // We need to set the options and check how many bins etc.
    seqan3::interleaved_bloom_filter<> ibf{seqan3::bin_count{bin_count},
                                           seqan3::bin_size{computed_bin_size},
                                           seqan3::hash_function_count{3}};

    size_t window_size = options.kmer_size - 1;
    for (size_t bin_id = 0; bin_id < bin_count; ++bin_id) {
        // make more efficient by providing a hasher.
        auto kmer_tree = forest[bin_id] | libjst::labelled()
                                        | libjst::coloured()
                                        | libjst::trim(window_size)
                                        | libjst::prune_unsupported()
                                        | libjst::left_extend(window_size)
                                        | libjst::merge();

        libjst::tree_traverser_base kmer_path{kmer_tree};
        for (auto it = kmer_path.begin(); it != kmer_path.end(); ++it) {
            auto label = *it;
            auto hash_seq = label.sequence() | seqan3::views::kmer_hash(seqan3::ungapped{options.kmer_size});
            for (uint64_t hash_value : hash_seq)
                ibf.emplace(hash_value, seqan3::bin_index{bin_id});
        }
    }

    return ibf;
}

} // namespace jstmap
