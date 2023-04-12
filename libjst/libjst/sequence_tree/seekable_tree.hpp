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
        using sink_type = libjst::tree_sink_t<base_tree_t>;
        using base_cargo_type = libjst::tree_label_t<base_tree_t>;
        using breakend_iterator = decltype(std::declval<base_node_type const &>().low_boundary().get_breakend());
        using difference_type = std::iter_difference_t<breakend_iterator>;

        class node_impl;
        class cargo_impl;

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
            base_node_type && base_root = libjst::root(_wrappee);
            difference_type root_index = std::ranges::distance(std::ranges::begin(data().variants()),
                                                               base_root.low_boundary().get_breakend());
            seek_position root_offset{};
            root_offset.reset(root_index, base_root.low_boundary().get_breakend_site());
            return node_impl{std::move(base_root), std::move(root_offset)};
        }

        constexpr sink_type sink() const noexcept {
            return libjst::sink(_wrappee);
        }

        constexpr auto const & data() const noexcept {
            return _wrappee.data();
        }

        constexpr node_impl seek(seek_position offset) const {
            return offset.visit([&] (auto path_descriptor) {
                breakend_iterator seek_breakend = std::ranges::next(std::ranges::begin(data().variants()),
                                                                    offset.get_variant_index());

                return unwind(std::move(path_descriptor), std::move(seek_breakend));
            });
        }

    private:

        constexpr node_impl unwind(breakpoint_end site, breakend_iterator seek_breakend) const {
            node_impl tmp = root();
            difference_type breakend_idx = std::ranges::distance(std::ranges::begin(data().variants()), seek_breakend);
            seek_position initial_position{};
            initial_position.reset(breakend_idx, site);
            tmp.reset(breakend_site<breakend_iterator>{std::move(seek_breakend), site}, std::move(initial_position));
            return tmp;
        }

        constexpr node_impl unwind(alternate_path_descriptor const & descriptor, breakend_iterator seek_breakend) const {
            node_impl tmp = root();
            --seek_breakend;
            difference_type breakend_idx = std::ranges::distance(std::ranges::begin(data().variants()), seek_breakend);
            breakpoint_end low_end = (*seek_breakend).get_breakpoint_end();
            seek_position initial_position{};
            initial_position.reset(breakend_idx, low_end);
            tmp.reset(breakend_site<breakend_iterator>{std::move(seek_breakend), low_end}, std::move(initial_position));

            for (auto it = std::ranges::begin(descriptor); it != std::ranges::end(descriptor); ++it) {
                if (*it) {
                    tmp = *tmp.next_alt();
                } else {
                    tmp = *tmp.next_ref();
                }
            }
            return tmp;
        }
    };

    template <typename base_tree_t>
    class seekable_tree_impl<base_tree_t>::node_impl : public base_node_type {
    private:

        friend seekable_tree_impl;

        seek_position _seek_offset{};

        explicit constexpr node_impl(base_node_type && base_node, seek_position seek_offset) noexcept :
            base_node_type{std::move(base_node)},
            _seek_offset{std::move(seek_offset)}
        {}

    public:

        constexpr node_impl() = default;

        constexpr base_node_type base() const & noexcept {
            return static_cast<base_node_type const &>(*this);
        }

        constexpr base_node_type base() && noexcept {
            return static_cast<base_node_type &&>(*this);
        }

        constexpr std::optional<node_impl> next_alt() const noexcept {
            return visit<true>(base_node_type::next_alt());
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            return visit<false>(base_node_type::next_ref());
        }

        constexpr cargo_impl operator*() const noexcept {
            return cargo_impl{this};
        }

    protected:

        template <typename breakend_site_t>
        constexpr void reset_low(breakend_site_t && low_boundary) {
            base_node_type::reset_low(std::move(low_boundary));
        }

        template <typename breakend_site_t>
        constexpr void reset(breakend_site_t && low_boundary, seek_position offset) {
            reset_offset(std::move(offset));
            reset_low(std::move(low_boundary));
        }

        constexpr void reset_offset(seek_position offset) noexcept {
            _seek_offset = std::move(offset);
        }

    private:

        constexpr difference_type variant_index() const noexcept {
            return _seek_offset.get_variant_index();
        }

        template <bool is_alt>
        constexpr std::optional<node_impl> visit(auto maybe_child) const {
            if (maybe_child) {
                seek_position child_offset{_seek_offset};
                if (this->on_alternate_path()) {
                    child_offset.next_alternate_node(is_alt);
                } else {
                    difference_type next_index = variant_index() + 1;
                    if constexpr (is_alt) {
                        child_offset.initiate_alternate_node(next_index);
                    } else {
                        child_offset.reset(next_index, maybe_child->low_boundary().get_breakend_site());
                    }
                }
                return node_impl{std::move(*maybe_child), std::move(child_offset)};
            } else {
                return std::nullopt;
            }
        }

        constexpr seek_position const & tell() const noexcept {
            return _seek_offset;
        }

        constexpr friend bool operator==(node_impl const & lhs, sink_type const & rhs) noexcept
        {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    template <typename base_tree_t>
    class seekable_tree_impl<base_tree_t>::cargo_impl : public base_cargo_type {
    private:

        friend seekable_tree_impl;

        node_impl const * _node{};

        explicit constexpr cargo_impl(node_impl const * node) noexcept :
            base_cargo_type{*static_cast<base_node_type const &>(*node)},
            _node{node}
        {}

    public:

        cargo_impl() = default;

        constexpr seek_position const & position() const noexcept {
            assert(_node != nullptr);
            return _node->tell();
        }
    };

    namespace _tree_adaptor {
        inline constexpr struct _seek
        {
            template <typename tree_t, typename ...args_t>
            constexpr auto operator()(tree_t && tree, args_t &&... args) const
                noexcept(std::is_nothrow_constructible_v<
                            seekable_tree_impl<std::remove_reference_t<tree_t>>, args_t...>)
                -> seekable_tree_impl<std::remove_reference_t<tree_t>, args_t...>
            {
                using adapted_tree_t = seekable_tree_impl<std::remove_reference_t<tree_t>, args_t...>;
                return adapted_tree_t{(tree_t &&)tree, (args_t &&)args...};
            }

            template <typename ...args_t>
            constexpr auto operator()(args_t &&... args) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, _seek, args_t...>)
                -> jst::contrib::closure_result_t<_seek, args_t...>
            { // we need to store the type that needs to be called later!
                return jst::contrib::make_closure(_seek{}, (args_t &&)args...);
            }
        } seek{};
    } // namespace _tree_adaptor

    using _tree_adaptor::seek;
}  // namespace libjst
