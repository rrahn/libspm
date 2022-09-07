// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the main entry point of the just_map creator.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <cereal/archives/binary.hpp>

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>

#include <jstmap/global/application_logger.hpp>
#include <jstmap/create/create_main.hpp>
// #include <jstmap/create/journaled_sequence_tree_builder.hpp>
#include <jstmap/create/load_sequence.hpp>
#include <jstmap/create/options.hpp>
// #include <jstmap/create/serialise_jst.hpp>
#include <jstmap/create/vcf_parser.hpp>

namespace jstmap
{

int create_main(seqan3::argument_parser & create_parser)
{
    create_options options{};

    create_parser.add_positional_option(options.sequence_file,
                                        "The input file.",
                                        seqan3::input_file_validator{{"fa", "fasta"}});
    create_parser.add_positional_option(options.output_file,
                                        "The output file.",
                                        seqan3::output_file_validator{seqan3::output_file_open_options::create_new,
                                                                      {"jst"}});
    create_parser.add_flag(options.is_quite,
                           '\0',
                           "quite",
                           "No logging output will be emitted.");
    create_parser.add_flag(options.is_verbose,
                           '\0',
                           "verbose",
                           "Verbose logging output will be emitted.");
    create_parser.add_option(options.vcf_file,
                             '\0',
                             "vcf",
                             "The vcf file to construct the index for. Note the path given to the sequence file "
                             "must contain the associated contigs for this vcf file.",
                             seqan3::option_spec::standard,
                             seqan3::input_file_validator{{"vcf"}});
    create_parser.add_option(options.bin_count,
                             'b',
                             "bin-count",
                             "The number of bins used in the partitioned jst.",
                             seqan3::option_spec::standard);

    try
    {
        create_parser.parse();
    }
    catch (seqan3::argument_parser_error const & ex)
    {
        get_application_logger()(verbosity_level::standard, logging_level::error,
                                 "While parsing command line options: ", ex.what());
        return -1;
    }

    // ----------------------------------------------------------------------------
    // Initialise the global logger
    // ----------------------------------------------------------------------------

    jstmap::application_logger logger{true, options.is_quite
                                                ? verbosity_level::quite
                                                : options.is_verbose
                                                    ? verbosity_level::verbose
                                                    : verbosity_level::standard};
    set_application_logger(&logger);
    jstmap::application_logger & log = get_application_logger();

    // ----------------------------------------------------------------------------
    // Run the index creation.
    // ----------------------------------------------------------------------------

    log(verbosity_level::standard, logging_level::info, "Start jst creation");

    int error_code = 0;
    try
    {
        if (create_parser.is_option_set("vcf")) // Construct from the vcf file.
        {
            log(verbosity_level::standard, logging_level::info,
                "Create from vcf ", options.vcf_file, " and contigs ", options.sequence_file);

            // auto jst_per_contig = construct_jst_from_vcf(options.sequence_file, options.vcf_file);
            construct_jst_from_vcf2(options.sequence_file, options.vcf_file, options.output_file);

            // log(verbosity_level::standard, logging_level::info,
            //     "Generated ", jst_per_contig.size(), " journal sequence tree(s).");

            // // Create jst per contig and give it a proper name.
            // size_t contig_idx{};
            // std::ranges::for_each(jst_per_contig, [&] (auto const & jst)
            // {
            //     std::filesystem::path filename_with_id = options.output_file;
            //     std::filesystem::path new_filename = filename_with_id.stem();
            //     new_filename += std::filesystem::path{std::to_string(contig_idx++)};
            //     new_filename += filename_with_id.extension();
            //     filename_with_id.replace_filename(new_filename);

            //     partitioned_jst_t partitioned_jst{std::addressof(jst), options.bin_count};

            //     log(verbosity_level::standard, logging_level::info, "Serialise jst ", filename_with_id);
            //     serialise(jst, partitioned_jst, options.output_file);
            // });
        }
        // else // Construct from the sequence alignment.
        // {
        //     log(verbosity_level::standard, logging_level::info, "Create by alignment <", options.sequence_file, ">");
        //     auto sequences = load_sequences(options.sequence_file);
        //     log(verbosity_level::verbose, logging_level::info, "Loaded ", sequences.size(), " sequences");
        //     auto [jst, partitioned_jst] = build_journaled_sequence_tree(std::move(sequences), options.bin_count);

        //     log(verbosity_level::verbose, logging_level::info, "Serialise jst <", options.output_file, ">");
        //     serialise(jst, partitioned_jst, options.output_file);
        // }
    }
    catch (std::exception const & ex)
    {
        log(verbosity_level::standard, logging_level::error, "While creating the jst: ", ex.what());
        error_code = -1;
    }
    log(verbosity_level::standard, logging_level::info, "Stop jst creation");
    return error_code;
}

} // namespace jstmap
