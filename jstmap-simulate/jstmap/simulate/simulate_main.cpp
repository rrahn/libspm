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
    simulate_parser.add_option(options.read_size,
                               's',
                               "read-size",
                               "The size of the reads.",
                               seqan3::option_spec::standard,
                               seqan3::arithmetic_range_validator{10u, 500u});

    simulate_parser.add_option(options.read_count,
                               'c',
                               "read-count",
                               "The number of reads to sample.",
                               seqan3::option_spec::standard,
                               seqan3::arithmetic_range_validator{100u, std::numeric_limits<uint32_t>::max()});

    simulate_parser.add_option(options.error_rate,
                               'e',
                               "error-rate",
                               "The relative error rate.",
                               seqan3::option_spec::standard,
                               seqan3::arithmetic_range_validator{0,1});

    try
    {
        simulate_parser.parse();
    }
    catch (seqan3::argument_parser_error const & ex)
    {
        std::cerr << "ERROR: " << ex.what() << "\n";
        return -1;
    }

    // Load the sequences.
    try
    {
        std::cout << "load the jst\n";
        auto rcs_store = load_jst(options.input_file);
        std::cout << "Sample reads\n";
        auto start = std::chrono::high_resolution_clock::now();

        read_sampler sampler{rcs_store};
        auto sampled_reads = sampler(options.read_count, options.read_size);

        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Sample time: "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
                  << " sec\n";

        std::cout << "Write out reads\n";
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
        std::cout << "Write out time: "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
                  << " sec\n";
    }
    catch (std::exception const & ex)
    {
        std::cerr << "ERROR: " << ex.what() << "\n";
        return -1;
    }

    return 0;
}

} // namespace jstmap
