// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides test mockup for testing the sequence trees.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <vector>
#include <set>
#include <string_view>

#include <libjst/coverage/bit_coverage.hpp>
#include <libjst/coverage/range_domain.hpp>
#include <libjst/utility/bit_vector.hpp>
#include <libjst/variant/breakpoint.hpp>
#include <libjst/variant/concept.hpp>

namespace jst::test
{

    template <typename position_t = libjst::breakpoint,
              typename insertion_t = std::string_view,
              typename deletion_t = int,
              typename coverage_t = libjst::bit_vector<>>
    struct variant {
        position_t position{};
        insertion_t insertion{};
        deletion_t deletion{};
        coverage_t coverage{};
        libjst::range_domain<uint32_t> domain{};

    private:
        friend position_t & tag_invoke(std::tag_t<libjst::position>, variant & me) {
            return me.position;
        }
        friend position_t tag_invoke(std::tag_t<libjst::position>, variant const & me) {
            return me.position;
        }
        friend position_t & tag_invoke(std::tag_t<libjst::left_breakpoint>, variant & me) {
            return libjst::position(me);
        }
        friend position_t tag_invoke(std::tag_t<libjst::left_breakpoint>, variant const & me) {
            return libjst::position(me);
        }
        friend insertion_t const & tag_invoke(std::tag_t<libjst::alt_sequence>, variant const & me) {
            return me.insertion;
        }
        friend insertion_t & tag_invoke(std::tag_t<libjst::alt_sequence>, variant & me) {
            return me.insertion;
        }
        friend deletion_t tag_invoke(std::tag_t<libjst::breakpoint_span>, variant const & me) {
            return me.deletion;
        }
        friend deletion_t & tag_invoke(std::tag_t<libjst::breakpoint_span>, variant & me) {
            return me.deletion;
        }
        friend coverage_t const & tag_invoke(std::tag_t<libjst::coverage>, variant const & me) {
            return me.coverage;
        }
        friend coverage_t & tag_invoke(std::tag_t<libjst::coverage>, variant & me) {
            return me.coverage;
        }

        constexpr friend std::weak_ordering operator<=>(variant const & lhs, variant const & rhs) noexcept {
            if (auto cmp = libjst::left_breakpoint(lhs) <=> libjst::left_breakpoint(rhs); cmp != 0) {
                return cmp;
            } else {
                return libjst::alt_kind(lhs) <=> libjst::alt_kind(rhs);
            }
        }
    };

    // static_assert(libjst::covered_sequence_variant<variant<>>);

    template <typename source_t>
    class mock_store {
        using variant_t = variant<libjst::breakpoint, source_t>;

    public:

        using source_type = source_t;
        using variant_map_type = std::multiset<variant_t>;

    private:
        source_type _source{};
        variant_map_type _map{};
        libjst::range_domain<uint32_t> _domain{};

    public:

        mock_store() = default;
        mock_store(source_type source, uint32_t const dom_size) : _source{std::move(source)}, _domain{0, dom_size}
        {
            libjst::bit_vector<> init{};
            init.resize(size(), true);
            _map.insert(variant_t{.position{0, libjst::breakpoint_end::left},
                                  .insertion{},
                                  .deletion{},
                                  .coverage{init},
                                  .domain{_domain}});
            uint32_t src_size = std::ranges::size(_source);
            _map.insert(variant_t{.position{src_size, libjst::breakpoint_end::right},
                                  .insertion{},
                                  .deletion{},
                                  .coverage{init},
                                  .domain{_domain}});
        }

        void insert(variant_t var) noexcept {
            assert(libjst::left_breakpoint(var).value() + libjst::breakpoint_span(var) <= std::ranges::ssize(source()));
            var.domain = _domain;
            assert(std::ranges::size(libjst::coverage(var)) <= size());
            _map.insert(std::move(var));
        }

        std::size_t size() const noexcept {
            return _domain.size();
        }

        source_type const & source() const noexcept {
            return _source;
        }

        variant_map_type const & variants() const noexcept {
            return _map;
        }
    };

}  // namespace jst::test
