// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the main entry point of the just_map viewer.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>

#include <jstmap/global/application_logger.hpp>
#include <jstmap/global/load_jst.hpp>
#include <jstmap/global/bam_writer.hpp>
#include <jstmap/linear/options.hpp>
#include <jstmap/global/tool_parser.hpp>

namespace jstmap
{

int linear_main(seqan3::argument_parser & linear_parser)
{
    linear_options options{};
    add_base_options(linear_parser, options);

    linear_parser.add_positional_option(options.rcsdb_file,
                                        "The rcsdb file containing the reference sequence information.",
                                        seqan3::input_file_validator{{"jst"}});
    linear_parser.add_positional_option(options.sam_file,
                                        "The sam file containing the mapping information of the reads aligned against the rcsdb.",
                                        seqan3::input_file_validator{{"sam", "bam"}});
    linear_parser.add_option(options.haplotype_index,
                             'H',
                             "haplotype",
                             "The index of the haplotype to extract the alignment file for.",
                             seqan3::option_spec::standard,
                             seqan3::arithmetic_range_validator<size_t>{0, std::numeric_limits<size_t>::max()});

    linear_parser.add_option(options.output_file,
                             'o',
                             "output",
                             "The file containing the linearised mapping information. "
                             "If not specified, the sam file is written to the standard output.",
                             seqan3::option_spec::standard,
                             seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create,
                                                           {"sam", "bam"}});

    try
    {
        linear_parser.parse();
        initialise_logging_level(options);
        log_info("Starting linearisation of sam file");
        // auto configure_rcsdb =
        //     just load rcsdb
        //     then construct rcsdb absolute position converter
        //     then construct rcsdb coverage checker
        // auto exec_pipeline =
        //     when_all(configure_rcsdb, just load sam file)
        //     then //IN: wrapped rcsdb + coverage checker
        //         auto chunked_sam_records = sam_file | views::chunk(cs)
        //         make_stream(chunked_sam_records)
        //         on_stream(get_scheduler(exec)) // fork each of the following calls are executed on tp context
        //         transform_stream chunked_sam_records | disable_uncovered(ht_idx) // check if alignment is covered by the haplotype: in variant subtree must have coverage, in reference must have coverage.
        //         transform_stream covered_chunked_sam_records | transform(convert_to_absolute_position(ht_idx))  // conversion into absolute position of given ht_idx for aligned reads
        //         for_each_stream write to synchronised bam file
        // run(exec_pipeline);
    } catch (seqan3::argument_parser_error const & ex) {
        log_err("Program terminates because of ", ex.what());
        return EXIT_FAILURE;
    }
    log_info("Successfully finished linearisation");
    return EXIT_SUCCESS;
}

} // namespace jstmap
