// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides functions to write the results.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <ranges>
#include <utility>

// #include <seqan3/core/detail/pack_algorithm.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/output.hpp>

#include <jstmap/search/write_results.hpp>

namespace jstmap
{

// struct sam_record
// {
//     std::string qname{"*"};
//     int16_t flag{};
//     std::string rname{"*"};
//     int32_t pos{0};
//     int32_t mapq{255};
//     std::string cigar{"*"};
//     std::string rnext{"*"};
//     int32_t pnext{0};
//     int64_t tlen{0};
//     std::string seq{"*"};
//     std::string qual{"*"};

//     auto as_record() const
//     {
//         using record_t = seqan3::record<seqan3::type_list<std::string, int32_t>,
//                                            seqan3::fields<field::id,   field::mapq>>;

//         return record_t{qname, mapq};

//         seqan3::detail::for_each([&] (seqan3::field field_id)
//         {
//             switch (field_id)
//             {
//                 case seqan3::field::id: seqan3::get<seqan3::field::id>(tmp) = qname; break;
//                 // case seqan3::field::flag: seqan3::get<seqan3::field::flag>(tmp) = flag; break;
//                 case seqan3::field::ref_id: seqan3::get<seqan3::field::ref_id>(tmp) = rname; break;
//                 case seqan3::field::ref_offset: seqan3::get<seqan3::field::ref_offset>(tmp) = pos; break;
//                 case seqan3::field::mapq: seqan3::get<seqan3::field::mapq>(tmp) = mapq; break;
//                 // case seqan3::field::alignment:
//                 // case seqan3::field::cigar: seqan3::get<seqan3::field::cigar>(tmp) = cigar; break;
//                 // case seqan3::field::mate: seqan3::get<seqan3::field::mate>(tmp) = std::tie(rnext, pnext, tlen); break;
//                 case seqan3::field::seq: seqan3::get<seqan3::field::seq>(tmp) = seq; break;
//                 case seqan3::field::qual: seqan3::get<seqan3::field::qual>(tmp) = qual; break;
//                 // case seqan3::field::tags:
//                 default: throw std::invalid_argument{"The given field is not supported"};
//             }
//         }, field_ids...);

//         return tmp{};
//     }
// };

// Process:
// a) query -> produce results for multiple reference genomes -> dump all data
//
//
void write_results(sam_file_t & sam_file,
                   std::vector<jstmap::search_match> const & matches,
                   seqan::StringSet<std::views::all_t<raw_sequence_t const &>> const & queries)
{
    auto alignment_pairs_view = matches | std::views::transform([&] (auto const & match)
    {
        return std::pair{match.sequence(), queries[match.query_id] | std::views::all};
    });

    auto align_cfg = seqan3::align_cfg::method_global{} |
                     seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}} |
                     seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10},
                                                        seqan3::align_cfg::extension_score{-1}} |
                     seqan3::align_cfg::output_sequence1_id{} |
                     seqan3::align_cfg::output_sequence2_id{} |
                     seqan3::align_cfg::output_alignment{} |
                     seqan3::align_cfg::output_begin_position{} |
                     seqan3::align_cfg::output_end_position{} |
                     seqan3::align_cfg::output_score{};

    for (auto && align_result : seqan3::align_pairwise(alignment_pairs_view, align_cfg)) // | seqan3::align_cfg::on_result{write_record});
    {
        assert(align_result.sequence1_id() < matches.size());
        sam_file.emplace_back(matches[align_result.sequence1_id()].hit_coordinate.position,
                              queries[matches[align_result.sequence2_id()].query_id],
                              align_result.alignment());
    }
}

}  // namespace jstmap
