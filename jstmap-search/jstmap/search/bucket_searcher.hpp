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
        std::vector<std::ptrdiff_t> _last_position{};
        double _error_rate;
        /*parameters*/
    public:

        bucket_searcher() = delete;
        explicit bucket_searcher(bucket_t bucket, double error_rate) noexcept :
            _bucket{std::move(bucket)},
            _error_rate{error_rate}
        {
            _last_position.resize(_bucket.needle_list.size(), -1);
        }

        template <typename callback_t>
        constexpr void operator()(callback_t && callback) {
            // instantiate pigeonhole filter.
            pigeonhole_filter filter{_bucket, _error_rate};
            // run pigeonhole filter on subtree
            // std::cout << "Run filter\n" << std::flush;
            filter([&] (auto && cargo, auto && finder, auto && needle_position) {
                if (position_available(cargo, finder, needle_position)) return;

                uint32_t seed_size = endPosition(finder) - beginPosition(finder);
                seed_verifier verifier{_bucket, _error_rate, seed_size};
                verifier(cargo, finder, needle_position, callback);
            });
        }

    private:

        template <typename cargo_t, typename finder_t, typename needle_hit_t>
        constexpr bool position_available(cargo_t const & cargo,
                                          finder_t const & finder,
                                          needle_hit_t const & hit) noexcept
        {
            std::ptrdiff_t label_offset = std::ranges::ssize(cargo.sequence()) - beginPosition(finder);
            std::ptrdiff_t global_begin_pos = std::ranges::ssize(cargo.path_sequence()) - label_offset - hit.offset;
            assert(global_begin_pos >= 0);

            if (_last_position[hit.index] != global_begin_pos) {
                _last_position[hit.index] = global_begin_pos;
                return false;
            }
            return true;
        }
    };

}  // namespace jstmap
