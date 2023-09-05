// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides right extender.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/utility/detail/multi_invocable.hpp>

#include <jstmap/global/match_position.hpp>
#include <jstmap/search/seed_prefix_extender.hpp>
#include <jstmap/search/seed_suffix_extender.hpp>

#include <libjst/sequence_tree/seek_position.hpp>

namespace jstmap
{
    template <typename bucket_t>
    class seed_verifier
    {
        bucket_t const & _bucket;
        double _error_rate{};
        size_t _seed_size{};
    public:
        seed_verifier(bucket_t const & bucket, double error_rate, size_t seed_size) noexcept :
            _bucket{bucket},
            _error_rate{error_rate},
            _seed_size{seed_size}
        {}

        template <typename cargo_t, typename finder_t, typename needle_hit_t,  typename callback_t>
        constexpr void operator()(cargo_t && cargo,
                                  finder_t && finder,
                                  needle_hit_t && needle_hit,
                                  [[maybe_unused]] callback_t && callback) const
        {
            auto const & needle = _bucket.needle_list[needle_hit.i1];
            uint32_t max_errors = get_error_count(needle);
            std::ranges::subrange needle_suffix{std::ranges::next(std::ranges::begin(needle), needle_hit.i2 + _seed_size),
                                                std::ranges::end(needle)};

            seed_suffix_extender right_extender{_bucket.base_tree, std::move(needle_suffix), max_errors};
            // what do we actually need?

            right_extender(cargo, finder, [&] ([[maybe_unused]] auto && right_cargo,
                                               [[maybe_unused]] auto && right_finder,
                                               [[maybe_unused]] int32_t right_errors) {
                assert(right_errors >= 0);
                assert(static_cast<uint32_t>(right_errors) <= max_errors);
                // we need to build the left extender!
                std::ranges::subrange needle_prefix{std::ranges::begin(needle),
                                                    std::ranges::next(std::ranges::begin(needle), needle_hit.i2)};
                seed_prefix_extender left_extender{_bucket.base_tree, std::move(needle_prefix), max_errors - right_errors};
                left_extender(cargo, finder, [&] ([[maybe_unused]] auto && left_cargo,
                                                  [[maybe_unused]] auto && left_finder,
                                                  [[maybe_unused]] int32_t total_errors){
                    // TODO: wrap in special position object!
                    std::ptrdiff_t left_start = to_forward_end(endPosition(left_finder));
                    libjst::seek_position left_seek = left_extender.to_forward_position(left_cargo.position());
                    libjst::seek_position merged_seek = merge(std::move(left_seek), right_cargo.position());

                    callback(needle_hit.i1,
                             match_position{.tree_position = std::move(merged_seek), .label_offset = left_start});
                });
            });

        }
    private:

        template <typename needle_t>
        constexpr uint32_t get_error_count(needle_t const & needle) const noexcept {
            return static_cast<uint32_t>(floor(_error_rate * length(needle)));
        }

        constexpr libjst::seek_position merge(libjst::seek_position left, libjst::seek_position right) const noexcept {
            return right.visit(seqan3::detail::multi_invocable{
                [&] (libjst::alternate_path_descriptor const & right_descriptor) {
                    for (bool elem : right_descriptor)
                        left.next_alternate_node(elem);

                    return left;
                },
                [&] (...) { return left; }
            });
        }

        constexpr std::ptrdiff_t to_forward_end(std::ptrdiff_t const reverse_position) const noexcept {
            return std::ranges::ssize(_bucket.base_tree.data().source()) - reverse_position;
        }
    };
}  // namespace jstmap
