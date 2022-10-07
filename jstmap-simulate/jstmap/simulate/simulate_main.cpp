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

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>
#include <seqan3/io/sequence_file/output.hpp>

#include <jstmap/create/serialise_jst.hpp>
#include <jstmap/global/jstmap_type_alias.hpp>
#include <jstmap/global/load_jst.hpp>
#include <jstmap/simulate/options.hpp>
#include <jstmap/simulate/simulate_main.hpp>

namespace jstmap
{

auto sample_reads(jst_t const & jst, size_t const read_count, size_t const read_size)
{
    using jst::contrib::operator""_dna4;

    auto enumerator = jst.context_enumerator(read_size);
    size_t const sample_rate = std::ceil(jst.reference_at(0).size() / read_count);
    size_t counter = 0;

    std::vector<raw_sequence_t> sampled_reads{};

    for (auto it = enumerator.begin(); it != enumerator.end(); ++it)
    {
        if (++counter % sample_rate == 0)
        {
            auto read = *it;
            if (std::ranges::any_of(read, [] (auto symbol) { return symbol == 'N'_dna4; }))
                continue;

            sampled_reads.emplace_back(read.begin(), read.end());
        }
    }

    return sampled_reads;
}

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
        auto jst = load_jst(options.input_file);
        std::cout << "Sample reads\n";
        auto sampled_reads = sample_reads(jst, options.read_count, options.read_size);
        seqan3::sequence_file_output sout{options.output_file};

        std::cout << "Write out reads\n";
        size_t sample_idx{};
        std::string sampled_read_id{"sample_read_"};
        std::ranges::for_each(sampled_reads, [&] (raw_sequence_t const & sequence)
        {
            sout.emplace_back(sequence, sampled_read_id + std::to_string(sample_idx++));
        });
    }
    catch (std::exception const & ex)
    {
        std::cerr << "ERROR: " << ex.what() << "\n";
        return -1;
    }

    return 0;
}

} // namespace jstmap
