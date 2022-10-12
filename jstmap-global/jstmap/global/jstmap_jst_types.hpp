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

#include <libcontrib/seqan/alphabet.hpp>

#include <libjst/journaled_sequence_tree/journaled_sequence_tree_forward.hpp>
#include <libjst/journaled_sequence_tree/journaled_sequence_tree_model.hpp>
#include <libjst/utility/bit_vector.hpp>
#include <libjst/sequence_variant/variant_snp.hpp>
#include <libjst/sequence_variant/variant_generic.hpp>
#include <libjst/sequence_variant/variant_store_composite.hpp>
#include <libjst/sequence_variant/variant_store_covered.hpp>

namespace jstmap
{

using alphabet_t = jst::contrib::dna5;
using snp_t = libjst::snp_variant<alphabet_t>;
using indel_t = libjst::generic_variant<alphabet_t>;
using coverage_t = libjst::bit_vector<>;
using base_sequence_t = std::vector<alphabet_t>;

using snp_store_t = libjst::variant_store_covered<std::vector<snp_t>, coverage_t>;
using indel_store_t = libjst::variant_store_covered<std::vector<indel_t>, coverage_t>;
using variant_store_t = libjst::variant_store_composite<snp_store_t, indel_store_t>;

using jst_model_t = libjst::journaled_sequence_tree_model<base_sequence_t, variant_store_t>;
using fwd_jst_t = libjst::journaled_sequence_tree_forward_<jst_model_t>;

}  // namespace jstmap
