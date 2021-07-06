// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides .
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <filesystem>

namespace jstmap
{

struct search_options
{
    std::filesystem::path jst_input_file_path{}; //!< The file path to the journaled sequence tree.
    std::filesystem::path query_input_file_path{}; //!< The file path containing the queries.
    std::filesystem::path index_input_file_path{}; //!< The file path containing the ibf index.
    std::filesystem::path map_output_file_path{}; //!< The file path to write the alignment map file to.
    float error_rate{}; //!< The error rate to use for mapping the reads.
    size_t thread_count{1}; //!< The number of threads to use for the program.
};

}  // namespace jstmap
