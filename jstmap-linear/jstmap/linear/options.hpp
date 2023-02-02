// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the options for the view sub-command.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <filesystem>

#include <jstmap/global/options_base.hpp>

namespace jstmap
{

struct linear_options : public options_base
{
    std::filesystem::path rcsdb_file{}; //!< The file path to the database containing the reference information.
    std::filesystem::path sam_file{}; //!< The path to the sam file containing the mapping information.
    std::filesystem::path output_file{}; //!< The path to the linearised output sam file.
    size_t haplotype_index{0}; //!< The haplotype index to linearise.
};

}  // namespace jstmap
