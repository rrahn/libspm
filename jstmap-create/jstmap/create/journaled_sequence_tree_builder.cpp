// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <iterator>
#include <tuple>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>

#include <jstmap/create/journaled_sequence_tree_builder.hpp>

namespace jstmap
{

template <typename reference_t, typename sequence_t>
auto compress(reference_t const & reference, sequence_t const & sequence)
{
    auto align_cfg = seqan3::align_cfg::method_global{} |
                     seqan3::align_cfg::scoring_scheme{jstmap::scoring_scheme{seqan3::match_score{5}, seqan3::mismatch_score{-4}}} |
                     seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10},
                                                        seqan3::align_cfg::extension_score{-1}} |
                     seqan3::align_cfg::output_sequence1_id{} |
                     seqan3::align_cfg::output_sequence2_id{} |
                     seqan3::align_cfg::output_alignment{} |
                     seqan3::align_cfg::output_begin_position{} |
                     seqan3::align_cfg::output_end_position{} |
                     seqan3::align_cfg::output_score{};

    auto align_range = seqan3::align_pairwise(std::tie(reference, sequence), align_cfg);
    auto result = *align_range.begin();
    return result.alignment();
}

std::pair<jst_t, partitioned_jst_t> build_journaled_sequence_tree(std::vector<raw_sequence_t> && sequences,
                                                                  const uint32_t bin_count /* = 1 */)
{
    assert(!sequences.empty());

    // Move ownership of first sequence to jst.
    jst_t jst{std::move(sequences[0])};

    // Add first sequence as empty alignment against themself.
    jst.add(compress(jst.reference_at(0), jst.reference_at(0)));

    // Align remaining sequences against reference sequence and add it to the jst.
    std::for_each(std::next(sequences.begin()), sequences.end(), [&] (auto const & sequence)
    {
        jst.add(compress(jst.reference_at(0), sequence));
    });

    // Build partitioned journaled sequence tree over the jst
    partitioned_jst_t partitioned_jst{std::addressof(jst), bin_count};

    return std::pair{std::move(jst), std::move(partitioned_jst)};
}

}  // namespace jstmap
