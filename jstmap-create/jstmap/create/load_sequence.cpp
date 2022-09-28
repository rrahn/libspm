// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/std/algorithm>

#include <seqan3/io/sequence_file/input.hpp>

#include <jstmap/create/load_sequence.hpp>

namespace jstmap
{

struct input_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = jst::contrib::dna5;
    using sequence_legal_alphabet = jst::contrib::dna15;
};

std::vector<raw_sequence_t> load_sequences(std::filesystem::path const & sequence_file)
{
    std::vector<raw_sequence_t> sequences{};

    seqan3::sequence_file_input<input_traits> fin{sequence_file.c_str()};
    std::ranges::for_each(fin, [&] (auto const & sequence_record)
    {
        sequences.push_back(seqan3::get<seqan3::field::seq>(sequence_record));
    });

    return sequences;
}

} // namespace jstmap
