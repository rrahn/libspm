// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the function for loading in a reference sequence for the simulated alignment.
 * \author Tom Lukas Lankenau <tom.lankenau AT fu-berlin.de>
 */

#include <stdexcept> // invalid_argument

#include <seqan3/io/sequence_file/input.hpp>
#include <filesystem>

#include <jstmap/simulate/load_reference.hpp>

namespace jstmap
{

raw_sequence_t load_reference(std::filesystem::path const & sequence_file)
{
    seqan3::sequence_file_input<sequence_input_traits> fin{sequence_file};
    auto it = fin.begin();
    if (it == fin.end())
        throw std::invalid_argument("Input file is empty.");

    return std::move((*it).sequence());
}

} // namespace jstmap
