// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/std/algorithm>

#include <seqan3/io/sequence_file/input.hpp>

#include <jstmap/index/load_sequence.hpp>

namespace jstmap
{

std::vector<raw_sequence_t> load_sequences(std::filesystem::path const & sequence_file)
{
    std::vector<raw_sequence_t> sequences{};

    seqan3::sequence_file_input fin{sequence_file.c_str()};
    std::ranges::for_each(fin, [&] (auto const & sequence_record)
    {
        sequences.push_back(seqan3::get<seqan3::field::seq>(sequence_record));
    });

    return sequences;
}

} // namespace jstmap
