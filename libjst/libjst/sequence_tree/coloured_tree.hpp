// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also availabel at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides journaled sequence tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <libcontrib/closure_adaptor.hpp>

#include <libjst/sequence_tree/concept.hpp>
#include <libjst/variant/concept.hpp>

namespace libjst
{

    template <typename wrapped_tree_t>
    class coloured_tree {
    private:
        using base_node_type = libjst::tree_node_t<wrapped_tree_t>;
        using sink_type = libjst::tree_sink_t<wrapped_tree_t>;
        using base_cargo_type = libjst::tree_label_t<wrapped_tree_t>;

        using boundary_type = typename base_node_type::position_type;
        using delta_reference = typename boundary_type::delta_reference;
        using coverage_type = libjst::variant_coverage_t<delta_reference>;

        class node_impl;
        class cargo_impl;

        wrapped_tree_t _wrappee{};
        coverage_type const * _coverage{};

    public:

        template <typename wrappee_t>
            requires (!std::same_as<std::remove_cvref_t<wrappee_t>, coloured_tree> &&
                      std::constructible_from<wrapped_tree_t, wrappee_t>)
        constexpr explicit coloured_tree(wrappee_t && wrappee) noexcept : _wrappee{(wrappee_t &&)wrappee}
        {
            _coverage = std::addressof(libjst::coverage(*(data().variants().begin())));
        }

        constexpr node_impl root() const noexcept {
            return node_impl{libjst::root(_wrappee), this};
        }

        constexpr sink_type sink() const noexcept {
            return libjst::sink(_wrappee);
        }

        constexpr auto const & data() const noexcept {
            return _wrappee.data();
        }
    };

    template <typename wrapped_tree_t>
    class coloured_tree<wrapped_tree_t>::node_impl : public base_node_type {
    private:

        friend coloured_tree;

        using base_t = base_node_type;

        coloured_tree const * _host{};

        explicit constexpr node_impl(base_t && base_node, coloured_tree const * host) noexcept :
                base_t{std::move(base_node)},
                _host{host}
        {}

    public:

        node_impl() = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            return visit_next(base_t::next_alt());
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            return visit_next(base_t::next_ref());
        }

        constexpr cargo_impl operator*() const noexcept {
            if (this->from_reference()) {
                return cargo_impl{*(static_cast<base_t const &>(*this)), _host->_coverage};
            } else {
                return cargo_impl{*(static_cast<base_t const &>(*this)),
                            std::addressof(libjst::coverage(*(this->low_boundary())))};
            }
        }

    private:

        constexpr std::optional<node_impl> visit_next(std::optional<base_t> && base_child) const noexcept {
            if (base_child.has_value())
                return node_impl{std::move(*base_child), _host};
            return std::nullopt;
        }

        constexpr friend bool operator==(node_impl const & lhs, sink_type const & rhs) noexcept
        {
            return static_cast<base_t const &>(lhs) == rhs;
        }
    };

    template <typename wrapped_tree_t>
    class coloured_tree<wrapped_tree_t>::cargo_impl : public base_cargo_type {
    private:
        [[no_unique_address]] coverage_type const * _coverage{};

        friend coloured_tree;

        constexpr explicit cargo_impl(base_cargo_type base_cargo, coverage_type const * coverage) :
            base_cargo_type{std::move(base_cargo)},
            _coverage{coverage}
        {}
    public:

        cargo_impl() = default;

        constexpr coverage_type const & coverage() const noexcept {
            return *_coverage;
        }
    };

    namespace _tree_adaptor {
        struct _coloured {

            template <typename tree_t, typename ...args_t>
            constexpr auto operator()(tree_t && tree, args_t &&... args) const
                noexcept(std::is_nothrow_constructible_v<coloured_tree<std::remove_reference_t<tree_t>>>)
                -> coloured_tree<std::remove_reference_t<tree_t>>
            {
                using adapted_tree_t = coloured_tree<std::remove_reference_t<tree_t>>;
                return adapted_tree_t{(tree_t &&)tree, (args_t &&)args...};
            }

            template <typename ...args_t>
            constexpr auto operator()(args_t &&... args) const
                noexcept(std::is_nothrow_invocable_v<libjst::tag_t<jst::contrib::make_closure>, _coloured, args_t...>)
                -> jst::contrib::closure_result_t<_coloured, args_t...>
            {
                return jst::contrib::make_closure(_coloured{}, (args_t &&)args...);
            }
        };

        inline constexpr _coloured coloured{};
    } // namespace _tree_adaptor

    using _tree_adaptor::coloured;
}  // namespace libjst
