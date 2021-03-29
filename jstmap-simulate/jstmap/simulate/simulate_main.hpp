// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <jstmap/simulate/global_types.hpp>

namespace seqan3
{

class argument_parser;

} // namespace seqan3

namespace jstmap
{

int simulate_main(seqan3::argument_parser &);

} // namespace jstmap
