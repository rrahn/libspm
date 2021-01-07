// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the main entry point of the just_map indexer.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <cereal/archives/binary.hpp>

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>

#include <jstmap/index/load_sequence.hpp>
#include <jstmap/index/index_main.hpp>
#include <jstmap/index/journaled_sequence_tree_builder.hpp>
#include <jstmap/index/options.hpp>

namespace jstmap
{

int index_main(seqan3::argument_parser & index_parser)
{
    index_options options{};

    index_parser.add_positional_option(options.input_file,
                                       "The input file.",
                                       seqan3::input_file_validator{{"fa", "fasta"}});
    index_parser.add_positional_option(options.output_file, "The output file.");

    try
    {
        index_parser.parse();
    }
    catch (seqan3::argument_parser_error const & ex)
    {
        std::cerr << "ERROR: " << ex.what() << "\n";
        return -1;
    }

    // Load the sequences.
    try
    {
        std::cout << "Loading sequences\n";
        auto sequences = load_sequences(options.input_file);

        libjst::journaled_sequence_tree tree = build_journaled_sequence_tree(std::move(sequences));

        std::ofstream output_stream{options.output_file.c_str()};
        cereal::BinaryOutputArchive binary_archive{output_stream};
        tree.save(binary_archive);
    }
    catch (std::exception const & ex)
    {
        std::cerr << "ERROR: " << ex.what() << "\n";
        return -1;
    }

    return 0;
}

} // namespace jstmap
