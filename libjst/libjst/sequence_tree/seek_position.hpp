// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seek position.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>

#include <cereal/types/base_class.hpp>

#include <seqan3/core/concept/cereal.hpp>

#include <libjst/sequence_tree/node_descriptor.hpp>
#include <libjst/sequence_tree/path_descriptor.hpp>

namespace libjst
{
    class seek_position {
    private:

        using index_t = uint_fast64_t;

        static constexpr uint32_t max_index_width_v{(sizeof(index_t) << 3) - 1};

        union {
            node_descriptor ref{};
            alternate_path_descriptor alt;
        } _descriptor{};

        index_t _active_descriptor : 1;
        index_t _variant_index : max_index_width_v;

    public:

        constexpr seek_position() noexcept : _active_descriptor{0}, _variant_index{0}
        {}

        constexpr void initiate_alternate_node(index_t const variant_index) noexcept {
            activate_alternate_node();
            _active_descriptor = 1;
            _variant_index = variant_index;
            alternate_node() = alternate_path_descriptor{};
        }

        constexpr void next_alternate_node(bool const is_alternate) noexcept {
            alternate_node().next();
            if (is_alternate) {
                alternate_node().set_alt();
            } else {
                alternate_node().set_ref();
            }
        }

        constexpr void reset(index_t const variant_index, node_descriptor const ref_state) noexcept {
            activate_reference_node();
            _variant_index = variant_index;
            reference_node() = ref_state;
        }

        constexpr index_t get_variant_index() const noexcept {
            return _variant_index;
        }

        template <typename visitor_t>
        constexpr auto visit(visitor_t && visitor) const {
            if (alternate_node_is_active()) {
                return std::invoke((visitor_t &&)visitor, alternate_node());
            } else {
                return std::invoke((visitor_t &&)visitor, reference_node());
            }
        }

        template <seqan3::cereal_output_archive archive_t>
        void save(archive_t & oarchive) const
        {
            oarchive(_active_descriptor, _variant_index);
            visit([&] (auto const & descriptor) { oarchive(descriptor); });
        }

        template <seqan3::cereal_input_archive archive_t>
        void load(archive_t & iarchive)
        {
            iarchive(_active_descriptor, _variant_index);
            visit([&] (auto const & descriptor) { iarchive(descriptor); });
        }

    protected:

        constexpr alternate_path_descriptor & alternate_node() noexcept {
            assert(alternate_node_is_active());
            return _descriptor.alt;
        }

        constexpr alternate_path_descriptor const & alternate_node() const noexcept {
            assert(alternate_node_is_active());
            return _descriptor.alt;
        }

        constexpr node_descriptor & reference_node() noexcept {
            assert(!alternate_node_is_active());
            return _descriptor.ref;
        }

        constexpr node_descriptor const & reference_node() const noexcept {
            assert(!alternate_node_is_active());
            return _descriptor.ref;
        }

        constexpr bool alternate_node_is_active() const noexcept {
            return _active_descriptor;
        }

        constexpr void activate_alternate_node() noexcept {
            _active_descriptor = 1;
        }

        constexpr void activate_reference_node() noexcept {
            _active_descriptor = 0;
        }
    };
}  // namespace libjst
