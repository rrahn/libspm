// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the serialiser function for the jst.
 * \author Tom Lukas Lankenau <tom.lankenau AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/filesystem>

#include <cereal/archives/binary.hpp>

namespace jstmap
{

template <typename rcs_store_t>
void serialise(rcs_store_t const & rcs_store, std::filesystem::path const & output_path)
{
    using namespace std::literals;

    std::ofstream output_stream{output_path.c_str()};
    if (!output_stream.good())
        throw std::runtime_error{"Couldn't open path for storing the rcs store! The path is ["s +
                                 output_path.string() +
                                 "]"s};

    cereal::BinaryOutputArchive binary_archive{output_stream};
    rcs_store.save(binary_archive);
}

} // namespace jstmap
