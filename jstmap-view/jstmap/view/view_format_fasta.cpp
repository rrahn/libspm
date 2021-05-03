// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <fstream>
#include <string>

#include <seqan3/io/sequence_file/output.hpp>

#include <jstmap/view/view_format_fasta.hpp>

namespace jstmap
{

void view_as_format(jst_t const & jst, size_t const haplotype_index)
{
    using namespace std::literals;

    seqan3::sequence_file_output view_out{std::cout, seqan3::format_fasta{}};
    view_out.emplace_back(jst.sequence_at(haplotype_index),
                          "ID_"s + std::to_string(haplotype_index));
}

} // namespace jstmap
