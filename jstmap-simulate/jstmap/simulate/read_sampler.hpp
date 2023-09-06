// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <random>
#include <tuple>

#include <libjst/sequence_tree/seek_position.hpp>
#include <libjst/sequence_tree/stats.hpp>

#include <jstmap/global/jstmap_types.hpp>
#include <jstmap/global/match_position.hpp>

namespace jstmap
{

    using read_type = reference_t;
    using sampled_read_type = std::tuple<read_type, match_position>;
    using sampled_read_list_type = std::vector<sampled_read_type>;
    class read_sampler {
    private:

        rcs_store_t const & _rcs_store;

    public:

        read_sampler(rcs_store_t const &);

        sampled_read_list_type operator()(size_t read_count, size_t read_size) const noexcept
        {
            libjst::tree_stats stats = compute_tree_stats(read_size);
            return sample_reads(stats, read_count, read_size);
        }

    private:
        sampled_read_list_type sample_reads(libjst::tree_stats const &, size_t const, size_t const) const;
        libjst::tree_stats compute_tree_stats(size_t const) const;
        std::vector<size_t> generate_sample_positions(size_t const, std::pair<size_t, size_t> const, std::mt19937 &) const;
        std::pair<size_t, size_t> compute_sample_range(libjst::tree_stats const &, size_t const) const noexcept;
        bool validate_sample(sampled_read_type const &) const noexcept;
    };

} // namespace jstmap
