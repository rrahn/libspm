// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the main entry point of the just_map simulator.
 * \author Tom Lukas Lankenau <tom.lankenau AT fu-berlin.de>
 */

#include <algorithm>
#include <chrono>

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>
#include <seqan3/io/sequence_file/output.hpp>

#include <jstmap/global/application_logger.hpp>
#include <jstmap/global/load_jst.hpp>
#include <jstmap/simulate/options.hpp>
#include <jstmap/simulate/simulate_main.hpp>
#include <jstmap/simulate/read_sampler.hpp>

namespace jstmap
{

int simulate_main(seqan3::argument_parser & simulate_parser)
{
    simulate_options options{};

    simulate_parser.add_positional_option(options.input_file,
                                          "The jst to sample reads from.",
                                          seqan3::input_file_validator{{"jst"}});
    simulate_parser.add_positional_option(options.output_file,
                                          "The file containing the sampled reads.",
                                          seqan3::output_file_validator{seqan3::output_file_open_options::create_new,
                                                                        {"fa", "fasta"}});

    simulate_parser.add_flag(options.is_quite,
                             'q',
                             "quite",
                             "Disables all logging.",
                             seqan3::option_spec::standard);
    simulate_parser.add_flag(options.is_verbose,
                             'v',
                             "verbose",
                             "Enables expansive debug logging.",
                             seqan3::option_spec::standard);

    simulate_parser.add_option(options.read_size,
                               's',
                               "read-size",
                               "The size of the reads.",
                               seqan3::option_spec::standard,
                               seqan3::arithmetic_range_validator{1u, 500u});

    simulate_parser.add_option(options.read_count,
                               'c',
                               "read-count",
                               "The number of reads to sample.",
                               seqan3::option_spec::standard,
                               seqan3::arithmetic_range_validator{1u, std::numeric_limits<uint32_t>::max()});

    simulate_parser.add_option(options.error_rate,
                               'e',
                               "error-rate",
                               "The relative error rate.",
                               seqan3::option_spec::standard,
                               seqan3::arithmetic_range_validator{0,1});

    try
    {
        simulate_parser.parse();
        if (options.is_quite) {
            get_application_logger().set_verbosity(verbosity_level::quite);
        } else if (options.is_verbose) {
            get_application_logger().set_verbosity(verbosity_level::verbose);
        }

        log_debug("Input file:", options.input_file.string());
        log_debug("Output file:", options.output_file.string());
        log_debug("Read size:", options.read_size);
        log_debug("Read count:", options.read_count);
        log_debug("Error rate:", options.error_rate);

    }
    catch (seqan3::argument_parser_error const & ex)
    {
        log(logging_level::error, ex.what());
        return -1;
    }

    // Load the sequences.
    auto global_start = std::chrono::high_resolution_clock::now();
    try
    {
        log_info("Starting simulation");
        log_debug("Load jst from file", options.input_file);
        auto rcs_store = load_jst(options.input_file);

        auto start = std::chrono::high_resolution_clock::now();
        log_debug("Initiate read sampler");
        read_sampler sampler{rcs_store};
        auto sampled_reads = sampler(options.read_count, options.read_size);

        auto end = std::chrono::high_resolution_clock::now();
        log_debug("Sampling time:", std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), "s");
        log_debug("Save sampled reads");
        start = std::chrono::high_resolution_clock::now();

        seqan3::sequence_file_output sout{options.output_file};
        size_t sample_idx{};
        std::ranges::for_each(sampled_reads, [&] (auto const & sampled_read)
        {
            using namespace std::literals;
            auto && [read_sequence, tree_position, label_offset] = sampled_read;
            std::string id_buffer = "jstsim|"s + std::to_string(sample_idx++) + "|"s + std::to_string(label_offset);
            sout.emplace_back(read_sequence, id_buffer);
        });

        end = std::chrono::high_resolution_clock::now();
        log_debug("Saving time:", std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), "s");
    }
    catch (std::exception const & ex)
    {
        log_err(ex.what());
        return -1;
    }
    auto global_end = std::chrono::high_resolution_clock::now();
    log_info("Fished simulation [", std::chrono::duration_cast<std::chrono::seconds>(global_end - global_start).count(), "s]");

    return 0;
}

} // namespace jstmap
