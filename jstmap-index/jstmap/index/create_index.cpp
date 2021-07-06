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

#include <seqan3/range/views/kmer_hash.hpp>

#include <libjst/journal_sequence_tree_context_enumerator.hpp>

#include <jstmap/index/create_index.hpp>

namespace jstmap
{

seqan3::interleaved_bloom_filter<> create_index(jst_t const & jst, index_options const & options)
{
    partitioned_jst_t pjst{std::addressof(jst), options.bin_size};
    size_t ibf_size = 2ull * 1024ull * 1024ull * 1024ull; // 2GiBi
    size_t computed_bin_size = ibf_size / (((pjst.bin_count() + 63) / 64) * 64);
    // We need to set the options and check how many bins etc.
    seqan3::interleaved_bloom_filter<> ibf{seqan3::bin_count{pjst.bin_count()},
                                           seqan3::bin_size{computed_bin_size},
                                           seqan3::hash_function_count{3}};

    for (size_t bin_id = 0; bin_id < pjst.bin_count(); ++bin_id)
    {
        libjst::detail::journal_sequence_tree_context_enumerator enumerator{pjst.bin_at(bin_id), options.kmer_size};
        for (auto && context : enumerator)
            for (uint64_t hash_value : context | seqan3::views::kmer_hash(seqan3::ungapped{options.kmer_size}))
                ibf.emplace(hash_value, seqan3::bin_index{bin_id});
    }

    return ibf;
}

} // namespace jstmap
