// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides search match aligner.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <seqan3/utility/views/slice.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>

#include <jstmap/search/match_aligner.hpp>

namespace jstmap
{
    search_match match_aligner::operator()(match_position pos) const
    {
        auto node = _reference_tree.seek(pos.tree_position);
        auto cargo = *node;
        auto ref_sequence = cargo.path_sequence();

        std::ptrdiff_t errors = std::ceil(_error_rate * std::ranges::ssize(_query_sequence));
        std::ptrdiff_t begin_position = std::max<std::ptrdiff_t>(0, pos.label_offset - errors);
        std::ptrdiff_t end_position = std::min(pos.label_offset + std::ranges::ssize(_query_sequence) + errors,
                                               std::ranges::ssize(ref_sequence));

        auto align_config = match_aligner::get_alignment_config(
            seqan3::match_score{4},
            seqan3::mismatch_score{-5},
            seqan3::align_cfg::open_score{-10},
            seqan3::align_cfg::extension_score{-1}
        );

        auto ref_segment = ref_sequence | seqan3::views::slice(begin_position, end_position);
        auto results = seqan3::align_pairwise(std::tie(_query_sequence, ref_segment), align_config);
        auto pairwise_align_result = *results.begin();
        return search_match{std::move(pos), alignment_result{std::move(pairwise_align_result)}};
    }

    match_aligner::ref_tree_type match_aligner::init(rcs_store_t const & rcs_store) {
        size_t window_size = std::ranges::size(_query_sequence) - 1;
        return rcs_store | libjst::make_volatile()
                         | libjst::labelled()
                        //  | libjst::coloured()
                         | libjst::trim(window_size)
                        //  | libjst::prune()
                         | libjst::left_extend(window_size)
                         | libjst::merge()
                         | libjst::seek();
    }

}  // namespace jstmap
