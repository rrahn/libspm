// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <ranges>
#include <string>
#include <type_traits>

namespace seqan
{
    inline void set(std::string & target, std::string const & source)
    {
       target = source;
    }

    template <typename value_t, typename allocator_t>
    inline void set(std::vector<value_t, allocator_t> & target, std::vector<value_t, allocator_t> const & source)
    {
        target.resize(source.size());
        std::ranges::copy(source, target.begin());
    }
} // namespace seqan
