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

#include <libjst/sequence_tree/seek_position.hpp>

#include <jstmap/global/application_logger.hpp>
#include <jstmap/global/match_position.hpp>
#include <jstmap/search/seed_prefix_extender.hpp>
#include <jstmap/search/seed_suffix_extender.hpp>


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
        constexpr void operator()(cargo_t && seed_cargo,
                                  finder_t && seed_finder,
                                  needle_hit_t && needle_hit,
                                  [[maybe_unused]] callback_t && callback) const
        {
            auto && needle = _bucket.needle_list[needle_hit.index];
            log_debug("Verify needle: ", needle_hit);
            log_debug("Seed position: ", seed_cargo.position());
            uint32_t max_errors = get_error_count(needle);
            std::ptrdiff_t suffix_start = needle_hit.offset + needle_hit.count; // Can this be larger than length of needle?
            std::ranges::subrange needle_suffix{std::ranges::next(std::ranges::begin(needle), suffix_start),
                                                std::ranges::end(needle)};
            seed_suffix_extender suffix_extender{_bucket.base_tree, std::move(needle_suffix), max_errors};
            suffix_extender(seed_cargo, seed_finder, [&] (match_position end_position,
                                                          [[maybe_unused]] int32_t suffix_errors) {
                assert(suffix_errors >= 0);
                assert(static_cast<uint32_t>(suffix_errors) <= max_errors);
                log_debug("Found valid suffix at: ", end_position);
                log_debug("Prepare prefix at seed: ", seed_cargo.position());

                std::ranges::subrange needle_prefix{std::ranges::begin(needle),
                                                    std::ranges::next(std::ranges::begin(needle), needle_hit.offset)};
                seed_prefix_extender prefix_extender{_bucket.base_tree, std::move(needle_prefix), max_errors - suffix_errors};
                prefix_extender(seed_cargo, seed_finder, [&] (match_position begin_position,
                                                              [[maybe_unused]] int32_t total_errors){
                    log_debug("Extend prefix at: ", seed_cargo.position());
                    begin_position.tree_position = join(begin_position.tree_position, end_position.tree_position);
                    callback(needle_hit.index, std::move(begin_position));
                });
            });
        }
    private:

        template <typename needle_t>
        constexpr uint32_t get_error_count(needle_t const & needle) const noexcept {
            return static_cast<uint32_t>(floor(_error_rate * length(needle)));
        }

        constexpr libjst::seek_position join(libjst::seek_position prefix_position,
                                             libjst::seek_position suffix_position) const noexcept {
            return suffix_position.visit(seqan3::detail::multi_invocable{
                [&] (libjst::alternate_path_descriptor const & suffix_position_descriptor) {
                    prefix_position.visit([&] <typename descriptor_t> (descriptor_t const &) {
                        if constexpr (std::same_as<descriptor_t, libjst::breakpoint_end>)
                            prefix_position.initiate_alternate_node(suffix_position.get_variant_index());
                    });

                    auto it = std::ranges::next(suffix_position_descriptor.begin(), 1);
                    for (; it != suffix_position_descriptor.end(); ++it)
                        prefix_position.next_alternate_node(*it);

                    return prefix_position;
                },
                [&] (...) { return prefix_position; }
            });
        }

        constexpr std::ptrdiff_t to_forward_end(std::ptrdiff_t const reverse_position) const noexcept {
            return std::ranges::ssize(_bucket.base_tree.data().source()) - reverse_position;
        }
    };
}  // namespace jstmap
