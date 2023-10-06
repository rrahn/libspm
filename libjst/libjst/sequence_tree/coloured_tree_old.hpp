// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides sequence tree with additional coverage checking.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <memory>

#include <libcontrib/closure_adaptor.hpp>

#include <libjst/sequence_tree/concept.hpp>
#include <libjst/sequence_tree/extendable_tree.hpp>
#include <libjst/sequence_tree/rcs_node_traits.hpp>

namespace libjst
{

    template <typename derived_t, typename extended_node_t>
    class node_coverage_extension {
    private:

        friend derived_t;

        using variant_type = typename rcs_node_traits<extended_node_t>::variant_type;
        using coverage_type = libjst::variant_coverage_t<variant_type>;

        // std::shared_ptr<coverage_type> _coverage{};
        coverage_type const * _coverage{};

        template <typename base_label_t>
        class label_impl;

    private:

        node_coverage_extension() = default;

        constexpr void initialise() { // TODO: replace automatic coverage generation.
            // what is the automatic coverage size?
            // _coverage = std::make_shared<coverage_type>(as_derived(*this).rcs_store().size(), true);
            _coverage = std::addressof(libjst::coverage(as_derived(*this).left_variant()));
        };

        constexpr node_coverage_extension notify(extended_node_t const &) const {
            return node_coverage_extension{*this}; // copy current extension into child.
        }

        static constexpr derived_t & as_derived(node_coverage_extension & me) noexcept {
            return static_cast<derived_t &>(me);
        }

        static constexpr derived_t const & as_derived(node_coverage_extension const & me) noexcept {
            return static_cast<derived_t const &>(me);
        }

    public:

        template <typename base_label_t>
        constexpr auto label(base_label_t && base_label) const noexcept {
            coverage_type const * cov_ptr{_coverage};
            if (as_derived(*this).is_alt_node()) {
                cov_ptr = std::addressof(libjst::coverage(as_derived(*this).left_variant())); //TODO
            }
            return label_impl<std::remove_cvref_t<base_label_t>>{(base_label_t &&) base_label, cov_ptr};
        }
    };

    template <typename derived_t, typename extended_node_t>
    template <typename base_label_t>
    class node_coverage_extension<derived_t, extended_node_t>::label_impl : public base_label_t {
    private:

        [[no_unique_address]] coverage_type const * _coverage{};

        friend node_coverage_extension;

        constexpr explicit label_impl(base_label_t base_label, coverage_type const * coverage) :
            base_label_t{std::move(base_label)},
            _coverage{coverage}
        {}

    public:

        label_impl() = default;

        constexpr coverage_type const & coverage() const noexcept {
            assert(_coverage != nullptr);
            return *_coverage;
        }
    };

    template <typename base_tree_t>
    using coloured_tree_impl = extendable_tree<base_tree_t, node_coverage_extension>;

    namespace _tree_adaptor {
        inline constexpr struct _coloured
        {
            template <typename tree_t, typename ...args_t>
            constexpr auto operator()(tree_t && tree, args_t &&... args) const
                noexcept(std::is_nothrow_constructible_v<
                            coloured_tree_impl<std::remove_reference_t<tree_t>>, args_t...>)
                -> coloured_tree_impl<std::remove_reference_t<tree_t>, args_t...>
            {
                using adapted_tree_t = coloured_tree_impl<std::remove_reference_t<tree_t>, args_t...>;
                return adapted_tree_t{(tree_t &&)tree, (args_t &&)args...};
            }

            template <typename ...args_t>
            constexpr auto operator()(args_t &&... args) const
                noexcept(std::is_nothrow_invocable_v<libjst::tag_t<jst::contrib::make_closure>, _coloured, args_t...>)
                -> jst::contrib::closure_result_t<_coloured, args_t...>
            { // we need to store the type that needs to be called later!
                return jst::contrib::make_closure(_coloured{}, (args_t &&)args...);
            }
        } coloured{};
    } // namespace _tree_adaptor

    using _tree_adaptor::coloured;
}  // namespace libjst
