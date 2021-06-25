// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a struct with the possible command line options.
 * \author Tom Lukas Lankenau <tom.lankenau AT fu-berlin.de>
 */

#pragma once

#include <filesystem>

namespace jstmap
{

struct simulate_options
{
    std::filesystem::path input_file{}; //!< The file path containing the sequences to simulate.
    std::filesystem::path output_file{}; //!< The file path to write the simulated jst to.
    size_t read_size = 100; //!< The sampled read size.
    size_t read_count = 1000; //!< The number of sampled reads.
    double error_rate = 0.0; //!< The relative rate of errors in the simulated sequence.
};

}  // namespace jstmap
