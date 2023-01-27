// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <utility>

#include <jstmap/global/search_query.hpp>
#include <jstmap/search/options.hpp>
#include <jstmap/search/type_alias.hpp>

namespace jstmap
{

std::pair<size_t, std::vector<bucket_type>>
filter_queries(std::vector<search_query> const &, search_options const &);

} // namespace jstmap
