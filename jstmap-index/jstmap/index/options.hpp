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

struct index_options
{
    std::filesystem::path jst_input_file{}; //!< The file path the jst to build the index for.
    std::filesystem::path output_file{}; //!< The file path to write the constructed index to.
    bool is_quite{false}; //!< Wether the index app should run in quite mode.
    bool is_verbose{false}; //!< Wether the index app should run in verbose mode.
    size_t bin_size = 10'000; //!< The size of a bin for the index construction.
    size_t bin_overlap = 500; //!< The size of the bin overlap for the ibf construction.
    uint8_t kmer_size = 25; //!< The kmer-size to use for the ibf creation.
};

}  // namespace jstmap
