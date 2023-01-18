// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides matching strategy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/matcher/horspool_matcher.hpp>
#include <libjst/traversal/state_oblivious_traverser.hpp>

#include <jstmap/search/match.hpp>
#include <jstmap/search/type_alias.hpp>

namespace jstmap
{

    class matching_operation {
    private:

    public:

        matching_operation() = default;

        template <typename callback_t>
        constexpr void operator()(haystack_type && haystack, bucket_type const & bucket, callback_t && callback) const {

            std::ranges::for_each(bucket, [&] (auto const & query) {
                auto && [idx, sequence] = query;
                libjst::horspool_matcher matcher{sequence};
                libjst::state_oblivious_traverser traverser{};
                traverser(haystack, matcher, [&] ([[maybe_unused]] auto const & finder, [[maybe_unused]] auto const & jst_context) {
                    callback(match{idx});
                });
            });
        }
    };

}  // namespace jstmap
