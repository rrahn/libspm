// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides globally defined type aliases for the jstmap tools.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <vector>

#include <cereal/types/vector.hpp>

#include <seqan3/io/sequence_file/input.hpp>

#include <libcontrib/seqan/alphabet.hpp>

#include <libjst/referentially_compressed_sequence_store/rcs_store.hpp>
#include <libjst/utility/bit_vector.hpp>
#include <libjst/variant/single_base_replacement_store.hpp>

namespace jstmap
{

using alphabet_t = jst::contrib::dna5;
using coverage_t = libjst::bit_vector<>;
using reference_t = std::vector<alphabet_t>;
using sequence_collection_t = std::vector<reference_t>;

using snv_store_t = libjst::single_base_replacement_store<alphabet_t>;
using rcs_store_t = libjst::rcs_store<reference_t, snv_store_t>;

using variant_store_t = snv_store_t;
using variant_t = std::ranges::range_value_t<variant_store_t>;

struct sequence_input_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = alphabet_t;
    using sequence_legal_alphabet = jst::contrib::dna15; // conversion?
};

using sequence_file_t = seqan3::sequence_file_input<sequence_input_traits>;
using sequence_record_t = std::ranges::range_value_t<sequence_file_t>;
using record_sequence_t = std::remove_cvref_t<decltype(std::declval<sequence_record_t const &>().sequence())>;

}  // namespace jstmap
