// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bucket searcher.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <jstmap/search/pigeonhole_filter.hpp>
#include <jstmap/search/seed_verifier.hpp>

namespace jstmap
{
    template <typename bucket_t>
    class bucket_searcher
    {
    private:

        bucket_t _bucket;
        double _error_rate;
        /*parameters*/
    public:

        bucket_searcher() = delete;
        explicit bucket_searcher(bucket_t bucket, double error_rate) noexcept :
            _bucket{std::move(bucket)},
            _error_rate{error_rate}
        {}

        template <typename callback_t>
        constexpr void operator()(callback_t && callback) const {
            // instantiate pigeonhole filter.
            pigeonhole_filter filter{_bucket, _error_rate};
            // run pigeonhole filter on subtree
            filter([&] (auto && cargo, auto && finder, auto && needle_position) {
                uint32_t seed_size = endPosition(finder) - beginPosition(finder);
                seed_verifier verifier{_bucket, _error_rate, seed_size};
                verifier(cargo, finder, needle_position, callback);
            });
        }
    };

}  // namespace jstmap
