// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides merge tree with maximal branch-free nodes.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/closure_adaptor.hpp>

#include <libjst/sequence_tree/concept.hpp>
#include <libjst/sequence_tree/left_extend_tree.hpp>

namespace libjst
{

    template <typename base_tree_t>
    class left_extend_tree_impl {
    private:

        using base_node_type = libjst::tree_node_t<base_tree_t>;
        using base_label_type = libjst::tree_label_t<base_tree_t>;

        class node_impl;
        class sink_impl;
        class label_impl;

        base_tree_t _wrappee{};
        size_t _offset{};

    public:

        template <typename wrappee_t, std::integral offset_t>
            requires (!std::same_as<std::remove_cvref_t<wrappee_t>, left_extend_tree_impl> &&
                      std::constructible_from<base_tree_t, wrappee_t>)
        constexpr explicit left_extend_tree_impl(wrappee_t && wrappee, offset_t offset) noexcept :
            _wrappee{(wrappee_t &&)wrappee},
            _offset{std::move(offset)}
        {}

        constexpr node_impl root() const noexcept {
            return node_impl{libjst::root(_wrappee), _offset};
        }
        constexpr sink_impl sink() const noexcept {
            return sink_impl{libjst::sink(_wrappee)};
        }
    };

    template <typename base_tree_t>
    class left_extend_tree_impl<base_tree_t>::node_impl : public base_node_type {
    private:

        friend left_extend_tree_impl;

        size_t _offset{};
        size_t _min_position{};

        explicit constexpr node_impl(base_node_type && base_node, size_t offset) noexcept :
            base_node_type{std::move(base_node)},
            _offset{std::move(offset)},
            _min_position{base_node_type::left_breakpoint().value()}
        {}

        explicit constexpr node_impl(base_node_type && base_node, size_t offset, size_t min_position) noexcept :
            base_node_type{std::move(base_node)},
            _offset{offset},
            _min_position{min_position}
        {}

    public:

        node_impl() = default;

        constexpr auto operator*() const noexcept {
            base_node_type::left_breakpoint();
            return label_impl{*static_cast<base_node_type const &>(*this), _offset, _min_position};
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
                return node_impl{std::move(*maybe_child), _offset, _min_position};
            } else {
                return std::nullopt;
            }
        }

        friend bool operator==(node_impl const & lhs, sink_impl const & rhs) noexcept {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    template <typename base_tree_t>
    class left_extend_tree_impl<base_tree_t>::sink_impl {
    private:
        friend left_extend_tree_impl;

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

    template <typename base_tree_t>
    class left_extend_tree_impl<base_tree_t>::label_impl : public base_label_type {
    private:

        friend left_extend_tree_impl;

        size_t _left_extension{};
        size_t _min_position{};

        constexpr explicit label_impl(base_label_type base_label, size_t left_extension, size_t min_position) noexcept :
            base_label_type{std::move(base_label)},
            _left_extension{left_extension},
            _min_position{min_position}
        {}

    public:

        using typename base_label_type::size_type;

        label_impl() = default;

        constexpr auto sequence(size_type left_pos, size_type right_pos = base_label_type::npos) const {
            assert(left_pos >= _min_position);
            assert(left_pos <= right_pos);
            left_pos = std::max<std::ptrdiff_t>(_min_position, left_pos - _left_extension);
            return base_label_type::sequence(left_pos, right_pos);
        }

        constexpr auto sequence() const {
            return sequence(_min_position);
        }
    };

    namespace _tree_adaptor {
        inline constexpr struct _left_extend
        {
            template <typename labelled_tree_t, std::integral left_extension_t>
            constexpr auto operator()(labelled_tree_t && tree, left_extension_t left_extension) const
                -> left_extend_tree_impl<labelled_tree_t>
            {
                assert(left_extension >= 0);
                return left_extend_tree_impl<labelled_tree_t>{(labelled_tree_t &&) tree, std::move(left_extension)};
            }

            template <std::integral left_extension_t>
            constexpr auto operator()(left_extension_t const left_extension) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, left_extension_t>)
                -> jst::contrib::closure_result_t<_left_extend, left_extension_t>
            { // we need to store the type that needs to be called later!
                assert(left_extension >= 0);
                return jst::contrib::make_closure(_left_extend{}, left_extension);
            }
        } left_extend{};
    } // namespace _tree_adaptor

    using _tree_adaptor::left_extend;
}  // namespace libjst
