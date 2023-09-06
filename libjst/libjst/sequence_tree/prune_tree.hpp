// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides prune tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/debug_stream.hpp>

#include <libcontrib/closure_adaptor.hpp>

#include <libjst/coverage/concept.hpp>

namespace libjst
{
    template <typename base_tree_t>
        // requires covered tree
    class prune_tree_impl {
    private:
        using base_node_type = libjst::tree_node_t<base_tree_t>;
        using sink_type = libjst::tree_sink_t<base_tree_t>;
        using base_cargo_type = libjst::tree_label_t<base_tree_t>;

        using boundary_type = typename base_node_type::position_type;
        using delta_reference = typename boundary_type::delta_reference;
        using coverage_type = libjst::variant_coverage_t<delta_reference>;

        class node_impl;
        class cargo_impl;

        base_tree_t _wrappee{};

    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr prune_tree_impl() = default; //!< Default.

        template <typename wrapped_tree_t>
            requires (!std::same_as<std::remove_cvref_t<wrapped_tree_t>, prune_tree_impl> &&
                      std::constructible_from<base_tree_t, wrapped_tree_t>)
        explicit constexpr prune_tree_impl(wrapped_tree_t && wrappee) noexcept :
            _wrappee{(wrapped_tree_t &&)wrappee}
        {}
        //!\}

        constexpr node_impl root() const noexcept {
            base_node_type base_root = libjst::root(_wrappee);
            auto base_coverage = (*base_root).coverage();
            return node_impl{std::move(base_root), std::move(base_coverage)};
        }

        constexpr sink_type sink() const noexcept {
            return libjst::sink(_wrappee);
        }

        constexpr auto const & data() const noexcept {
            return _wrappee.data();
        }
   };

    template <typename base_tree_t>
    class prune_tree_impl<base_tree_t>::node_impl : public base_node_type {
    private:

        friend prune_tree_impl;

        coverage_type _path_coverage{};

        explicit constexpr node_impl(base_node_type && base_node, coverage_type coverage) noexcept :
            base_node_type{std::move(base_node)},
            _path_coverage{std::move(coverage)}
        {}

    public:

        node_impl() = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            return visit<true>(base_node_type::next_alt());
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            return visit<false>(base_node_type::next_ref());
        }

        constexpr cargo_impl operator*() const noexcept {
            return cargo_impl{*(static_cast<base_node_type const &>(*this)), std::addressof(_path_coverage)};
        }

    private:

        template <bool is_alt, typename maybe_child_t>
        constexpr std::optional<node_impl> visit(maybe_child_t maybe_child) const {
            if (maybe_child) {
                if (auto new_cov = compute_child_coverage<is_alt>(*maybe_child); new_cov.any()) {
                    node_impl new_child{std::move(*maybe_child), std::move(new_cov)};
                    return new_child;
                }
            }
            return std::nullopt;
        }

        template <bool is_alt>
        constexpr coverage_type compute_child_coverage(base_node_type const & base_child) const {
            if constexpr (is_alt) {
                // seqan3::debug_stream << "Cld cov: " << (*base_child).coverage() << "\n";
                // new_coverage &= (*base_child).coverage();
                return libjst::coverage_intersection(_path_coverage, (*base_child).coverage());
            } else if (this->on_alternate_path() && this->high_boundary().is_low_end()) {
                    // seqan3::debug_stream << "Alt cov:" << libjst::coverage(base_node_type::right_variant()) << "\n";
                    return libjst::coverage_difference(_path_coverage, libjst::coverage(*(this->high_boundary())));
            } else {
                return _path_coverage;
            }
            // seqan3::debug_stream << "New cov: " << new_coverage << "\n";
            // return new_coverage;
        }

        constexpr friend bool operator==(node_impl const & lhs, sink_type const & rhs) noexcept
        {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    template <typename base_tree_t>
    class prune_tree_impl<base_tree_t>::cargo_impl : public base_cargo_type {
    private:
        [[no_unique_address]] coverage_type const * _path_coverage{};

        friend prune_tree_impl;

        constexpr explicit cargo_impl(base_cargo_type base_cargo, coverage_type const * path_coverage) :
            base_cargo_type{std::move(base_cargo)},
            _path_coverage{path_coverage}
        {}
    public:

        cargo_impl() = default;

        constexpr coverage_type const & coverage() const noexcept {
            return *_path_coverage;
        }
    };


    namespace _tree_adaptor {
        inline constexpr struct _prune
        {
            template <typename covered_tree_t, typename ...args_t>
            constexpr auto operator()(covered_tree_t && tree, args_t &&... args) const
                noexcept(std::is_nothrow_constructible_v<
                            prune_tree_impl<std::remove_reference_t<covered_tree_t>>, args_t...>)
                -> prune_tree_impl<std::remove_reference_t<covered_tree_t>, args_t...>
            {
                using adapted_tree_t = prune_tree_impl<std::remove_reference_t<covered_tree_t>, args_t...>;
                return adapted_tree_t{(covered_tree_t &&)tree, (args_t &&)args...};
            }

            template <typename ...args_t>
            constexpr auto operator()(args_t &&... args) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, args_t...>)
                -> jst::contrib::closure_result_t<_prune, args_t...>
            { // we need to store the type that needs to be called later!
                return jst::contrib::make_closure(_prune{}, (args_t &&)args...);
            }
        } prune{};
    } // namespace _tree_adaptor

    using _tree_adaptor::prune;
}  // namespace libjst
