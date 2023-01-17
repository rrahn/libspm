// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seekable tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/closure_adaptor.hpp>

#include <libjst/sequence_tree/concept.hpp>
#include <libjst/sequence_tree/node_descriptor.hpp>
#include <libjst/sequence_tree/seek_position.hpp>

namespace libjst
{
    template <typename base_tree_t>
    class seekable_tree_impl {
    private:

        using base_node_type = libjst::tree_node_t<base_tree_t>;

        class node_impl;
        class sink_impl;
        class label_impl;

        base_tree_t _wrappee{};

    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr seekable_tree_impl() = default; //!< Default.

        template <typename wrapped_tree_t>
            requires (!std::same_as<wrapped_tree_t, seekable_tree_impl> &&
                      std::constructible_from<base_tree_t, wrapped_tree_t>)
        explicit constexpr seekable_tree_impl(wrapped_tree_t && wrappee) noexcept :
            _wrappee{(wrapped_tree_t &&)wrappee}
        {}
        //!\}

        constexpr node_impl root() const noexcept {
            return node_impl{libjst::root(_wrappee)};
        }

        constexpr sink_impl sink() const noexcept {
            return sink_impl{libjst::sink(_wrappee)};
        }

        constexpr node_impl seek(seek_position offset) const noexcept {
            auto node = root();
            node.seek(std::move(offset));
            return node;
        }
    };

    template <typename base_tree_t>
    class seekable_tree_impl<base_tree_t>::node_impl : public base_node_type {
    private:

        friend seekable_tree_impl;

        seek_position _seek_offset{};

        explicit constexpr node_impl(base_node_type && base_node) noexcept :
            base_node_type{std::move(base_node)}
        {
            _seek_offset.reset(variant_index(), static_cast<node_descriptor const &>(*this));
        }

        explicit constexpr node_impl(base_node_type && base_node, seek_position seek_pos) noexcept :
            base_node_type{std::move(base_node)},
            _seek_offset{std::move(seek_pos)}
        {}

    public:

        constexpr node_impl() = default;
        constexpr node_impl(node_impl const &) = default;
        constexpr node_impl(node_impl &&) = default;
        constexpr node_impl & operator=(node_impl const &) = default;
        constexpr node_impl & operator=(node_impl &&) = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            return visit<true>(base_node_type::next_alt());
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            return visit<false>(base_node_type::next_ref());
        }

        constexpr auto operator*() const noexcept {
            return label_impl{*static_cast<base_node_type const &>(*this), _seek_offset};
        }

    private:

        template <bool is_alt>
        constexpr std::optional<node_impl> visit(auto maybe_child) const {
            if (maybe_child) {
                node_impl child{std::move(*maybe_child), _seek_offset};
                if (base_node_type::on_alternate_path()) {
                    child._seek_offset.next_alternate_node(is_alt);
                } else {
                    if constexpr (is_alt) {
                       child._seek_offset.initiate_alternate_node(child.variant_index());
                    } else {
                        child._seek_offset.reset(child.variant_index(), static_cast<node_descriptor const &>(*this));
                    }
                }
                return child;
            } else {
                return std::nullopt;
            }
        }

        constexpr void seek(seek_position offset) {
            auto root_it = std::ranges::begin(base_node_type::rcs_store().variants());
            auto var_it = std::ranges::next(root_it, offset.get_variant_index());

            base_node_type::set_left(var_it);
            base_node_type::set_right(var_it);
            base_node_type::set_next(base_node_type::next_variant_after(base_node_type::get_right()));
            _seek_offset = std::move(offset);
            *this = offset.visit([&] (auto const & descriptor) { return unwind(descriptor); });
        }

        constexpr node_impl unwind(node_descriptor const & descriptor) {
            base_node_type::activate_state(static_cast<node_state>(descriptor));
            return *next_ref();
        }

        constexpr node_impl unwind(alternate_path_descriptor const & descriptor) {
            base_node_type::activate_state(node_state::branching_after_left_end);
            node_impl tmp{*next_alt()};

            auto it = std::ranges::next(std::ranges::begin(descriptor));
            for (; it != std::ranges::end(descriptor); ++it) {
                if (*it) {
                    tmp = *tmp.next_alt();
                } else {
                    tmp = *tmp.next_ref();
                }
            }
            return tmp;
        }

        constexpr seek_position tell() const noexcept {
            return _seek_offset;
        }

        constexpr auto variant_index() const noexcept {
            return std::ranges::distance(std::ranges::begin(base_node_type::rcs_store().variants()),
                                         base_node_type::get_left());
        }

        constexpr bool is_root() const noexcept {
            return base_node_type::get_left() == std::ranges::begin(base_node_type::rcs_store().variants());
        }
        constexpr friend bool operator==(node_impl const & lhs, sink_impl const & rhs) noexcept {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    template <typename base_tree_t>
    class seekable_tree_impl<base_tree_t>::label_impl : public libjst::tree_label_t<base_tree_t> {
    private:

        using base_label_t = libjst::tree_label_t<base_tree_t>;

        friend seekable_tree_impl;

        seek_position _seek_position{};

        template <typename base_label_t>
        explicit constexpr label_impl(base_label_t base_label, seek_position seek_position) noexcept :
            base_label_t{std::move(base_label)},
            _seek_position{std::move(seek_position)}
        {}

    public:

        label_impl() = default;

        constexpr seek_position const & position() const noexcept {

            return _seek_position;
        }
    };

    template <typename base_tree_t>
    class seekable_tree_impl<base_tree_t>::sink_impl {
    private:
        friend seekable_tree_impl;

        using base_sink_type = libjst::tree_sink_t<base_tree_t>;
        base_sink_type _base_sink{};

        constexpr explicit sink_impl(base_sink_type base_sink) : _base_sink{std::move(base_sink)}
        {}

        friend bool operator==(sink_impl const & lhs, base_node_type const & rhs) noexcept {
            return lhs._base_sink == rhs;
        }

    public:
        sink_impl() = default;
    };

    namespace _tree_adaptor {
        inline constexpr struct _seek
        {
            template <typename covered_tree_t, typename ...args_t>
            constexpr auto operator()(covered_tree_t && tree, args_t &&... args) const
                noexcept(std::is_nothrow_constructible_v<
                            seekable_tree_impl<std::remove_reference_t<covered_tree_t>>, args_t...>)
                -> seekable_tree_impl<std::remove_reference_t<covered_tree_t>, args_t...>
            {
                using adapted_tree_t = seekable_tree_impl<std::remove_reference_t<covered_tree_t>, args_t...>;
                return adapted_tree_t{(covered_tree_t &&)tree, (args_t &&)args...};
            }

            template <typename ...args_t>
            constexpr auto operator()(args_t &&... args) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, args_t...>)
                -> jst::contrib::closure_result_t<_seek, args_t...>
            { // we need to store the type that needs to be called later!
                return jst::contrib::make_closure(_seek{}, (args_t &&)args...);
            }
        } seek{};
    } // namespace _tree_adaptor

    using _tree_adaptor::seek;
}  // namespace libjst
