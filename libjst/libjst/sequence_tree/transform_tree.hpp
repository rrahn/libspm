// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides generalised transform tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>

#include <libcontrib/copyable_box.hpp>

#include <libjst/sequence_tree/concept.hpp>

namespace libjst
{
    template <typename base_tree_t, typename fn_t>
    class transform_tree_impl {
    private:

        using base_node_type = libjst::tree_node_t<base_tree_t>;

        class node_impl;
        class sink_impl;

        using tree_box_t = jst::contrib::copyable_box<base_tree_t>;
        using transform_box_t = jst::contrib::copyable_box<fn_t>;

        base_tree_t _wrappee{};
        transform_box_t _label_fn{};

    public:

        template <typename wrappee_t, typename label_transform_fn_t>
            requires (!std::same_as<std::remove_cvref_t<wrappee_t>, transform_tree_impl> &&
                      std::constructible_from<base_tree_t, wrappee_t>)
        constexpr explicit transform_tree_impl(wrappee_t && wrappee, label_transform_fn_t && fn) noexcept :
            _wrappee{(wrappee_t &&)wrappee},
            _label_fn{(label_transform_fn_t &&) fn}
        {}

        constexpr node_impl root() const noexcept {
            return node_impl{libjst::root(_wrappee), _label_fn};
        }
        constexpr sink_impl sink() const noexcept {
            return sink_impl{libjst::sink(_wrappee)};
        }
    };

    template <typename base_tree_t, typename fn_t>
    class transform_tree_impl<base_tree_t, fn_t>::node_impl : public base_node_type {
    private:

        friend transform_tree_impl;

        transform_box_t _fn{};

        explicit constexpr node_impl(base_node_type && base_node, transform_box_t fn) noexcept :
            base_node_type{std::move(base_node)},
            _fn{std::move(fn)}
        {}

    public:

        node_impl() = default;
        node_impl(node_impl const &) = default;
        node_impl(node_impl &&) = default;
        node_impl & operator=(node_impl const &) = default;
        node_impl & operator=(node_impl &&) = default;

        constexpr auto operator*() const noexcept {
            return std::invoke(_fn, *static_cast<base_node_type const &>(*this));
        }

        constexpr std::optional<node_impl> next_alt() const noexcept {
            return visit(base_node_type::next_alt());
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            return visit(base_node_type::next_ref());
        }

    private:

        constexpr std::optional<node_impl> visit(auto maybe_child) const {
            if (maybe_child) {
                return node_impl{std::move(*maybe_child), _fn};
            } else {
                return std::nullopt;
            }
        }

        friend bool operator==(node_impl const & lhs, sink_impl const & rhs) noexcept {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    template <typename base_tree_t, typename fn_t>
    class transform_tree_impl<base_tree_t, fn_t>::sink_impl {
    private:
        friend transform_tree_impl;

        using base_sink_type = libjst::tree_sink_t<base_tree_t>;
        base_sink_type _base_sink{};

        constexpr explicit sink_impl(base_sink_type base_sink) : _base_sink{std::move(base_sink)}
        {}

    private:

        sink_impl() = default;

        friend bool operator==(sink_impl const & lhs, base_node_type const & rhs) noexcept {
            return lhs._base_sink == rhs;
        }

    };

    namespace _tree_adaptor {
        inline constexpr struct _transform
        {
            template <typename covered_tree_t, typename fn_t>
            constexpr auto operator()(covered_tree_t && tree, fn_t && transform_fn) const
                noexcept(std::is_nothrow_constructible_v<
                            transform_tree_impl<std::remove_reference_t<covered_tree_t>, std::remove_reference_t<fn_t>>>)
                -> transform_tree_impl<std::remove_reference_t<covered_tree_t>, std::remove_reference_t<fn_t>>
            {
                using adapted_tree_t = transform_tree_impl<std::remove_reference_t<covered_tree_t>, std::remove_reference_t<fn_t>>;
                return adapted_tree_t{(covered_tree_t &&)tree, (fn_t &&) transform_fn};
            }

            template <typename fn_t>
            constexpr auto operator()(fn_t && transform_fn) const
                noexcept(std::is_nothrow_invocable_v<libjst::tag_t<jst::contrib::make_closure>, fn_t>)
                -> jst::contrib::closure_result_t<_transform, fn_t>
            { // we need to store the type that needs to be called later!
                return jst::contrib::make_closure(_transform{}, (fn_t &&)transform_fn);
            }
        } transform{};
    } // namespace _tree_adaptor

    using _tree_adaptor::transform;
}  // namespace libjst
