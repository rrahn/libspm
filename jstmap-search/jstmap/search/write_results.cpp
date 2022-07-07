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

#include <atomic>
#include <future>
#include <ranges>
#include <utility>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/contrib/parallel/buffer_queue.hpp>

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
void write_results(std::vector<search_match> const & matches,
                   std::vector<raw_sequence_t> const & queries,
                   search_options const & options)
{
    // ----------------------------------------------------------------------------
    // Configure the output record for the sam file.
    // ----------------------------------------------------------------------------

    using field_types_t = seqan3::type_list<int32_t,
                                            std::views::all_t<raw_sequence_t const &>,
                                            std::vector<seqan3::cigar>>;
    using field_ids_t = seqan3::fields<seqan3::field::ref_offset,
                                       seqan3::field::seq,
                                       seqan3::field::cigar>;
    using record_t  = seqan3::record<field_types_t, field_ids_t>;

    // ----------------------------------------------------------------------------
    // Configure the concurrent resources.
    // ----------------------------------------------------------------------------

    seqan3::contrib::dynamic_buffer_queue<record_t> result_queue{};

    // One consumer and thread_count - 1 many producer.
    bool const is_single_threaded = options.thread_count == 1;
    uint32_t const producer_count = (is_single_threaded) ? 1 : options.thread_count - 1;
    std::atomic<size_t> alignment_counter{}; // All alignments need to be computed.

    // ----------------------------------------------------------------------------
    // Define the producer job.
    // ----------------------------------------------------------------------------

    // FIXME: What if one alignment fails and does not produce a result? Then we never finish!
    auto async_push = [&] (auto const & align_result)
    {
        std::vector cigar = seqan3::detail::get_cigar_vector(align_result.alignment(),
                                                             align_result.sequence2_begin_position(),
                                                             align_result.sequence2_end_position());

        record_t record{matches[align_result.sequence1_id()].hit_coordinate.position,
                        queries[matches[align_result.sequence2_id()].query_id],
                        std::move(cigar)};
        result_queue.wait_push(std::move(record));
        ++alignment_counter;
    };

    // ----------------------------------------------------------------------------
    // Define the consumer job.
    // ----------------------------------------------------------------------------

    // If we have only one thread available, we have to first produce all alingments and then
    // write them to the file.
    std::launch launch_policy = (is_single_threaded) ? std::launch::deferred : std::launch::async;

    // Start thread to asynchronously write results into global sam file resource.
    std::future<void> serialiser = std::async(launch_policy, [&] ()
    {
        seqan3::sam_file_output sam_file{options.map_output_file_path, field_ids_t{}};

        while (true)
        {
            record_t record{};
            if (auto status = result_queue.wait_pop(record); status == seqan3::contrib::queue_op_status::closed)
                break;

            sam_file.push_back(std::move(record));
        }

        assert(result_queue.is_empty());
    });

    // ----------------------------------------------------------------------------
    // Configure and run the alignment.
    // ----------------------------------------------------------------------------

    auto align_cfg = seqan3::align_cfg::method_global{} |
                     seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}} |
                     seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10},
                                                        seqan3::align_cfg::extension_score{-1}} |
                     seqan3::align_cfg::output_sequence1_id{} |
                     seqan3::align_cfg::output_sequence2_id{} |
                     seqan3::align_cfg::output_alignment{} |
                     seqan3::align_cfg::output_begin_position{} |
                     seqan3::align_cfg::output_end_position{} |
                     seqan3::align_cfg::output_score{} |
                     seqan3::align_cfg::parallel{producer_count} |
                     seqan3::align_cfg::on_result{async_push};

    auto alignment_pairs_view = matches | std::views::transform([&] (auto const & match)
    {
        return std::pair{match.sequence(), queries[match.query_id] | std::views::all};
    });

    seqan3::align_pairwise(alignment_pairs_view, align_cfg);

    // ----------------------------------------------------------------------------
    // Closing step.
    // ----------------------------------------------------------------------------

    // Wait for the last produced alignment.
    while (alignment_counter.load() < matches.size())
    {}

    // Close the queue to signal that there are no more alignments coming.
    assert(alignment_counter.load() == matches.size());
    result_queue.close();

    // Wait for the async serialiser to write the last records.
    serialiser.get();
}

}  // namespace jstmap
