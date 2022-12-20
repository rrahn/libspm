// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides trim tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/closure_adaptor.hpp>
#include <libcontrib/copyable_box.hpp>

#include <libjst/sequence_tree/concept.hpp>
namespace libjst
{
    template <typename base_tree_t>
    class trim_tree_impl {
    private:
        using wrappee_t = jst::contrib::copyable_box<base_tree_t>;
        using base_node_type = libjst::tree_node_t<base_tree_t>;

        class node_impl;
        class sink_impl;

        wrappee_t _wrappee{};
        std::size_t _max_branch_size{};

    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr trim_tree_impl() = default; //!< Default.

        template <typename wrapped_tree_t>
        explicit constexpr trim_tree_impl(wrapped_tree_t && wrappee, std::size_t const max_branch_size) noexcept :
            _wrappee{(wrapped_tree_t &&)wrappee},
            _max_branch_size{max_branch_size}
        {}
        //!\}

        constexpr node_impl root() const noexcept {
            return node_impl{libjst::root(*_wrappee), _max_branch_size};
        }

        constexpr sink_impl sink() const noexcept {
            return sink_impl{libjst::sink(*_wrappee)};
        }
   };

    template <typename base_tree_t>
    class trim_tree_impl<base_tree_t>::node_impl : public base_node_type {
    private:

        friend trim_tree_impl;

        using difference_type = std::ptrdiff_t;

        std::size_t _max_branch_size{};
        difference_type _remaining_branch_size{};
        // bool _on_variant_branch{false};

        explicit constexpr node_impl(base_node_type base_node, std::size_t max_branch_size) noexcept :
            node_impl{std::move(base_node), max_branch_size, static_cast<difference_type>(max_branch_size)}
        {}

        explicit constexpr node_impl(base_node_type base_node,
                                     std::size_t max_branch_size,
                                     difference_type remaining_branch_size) noexcept :
            base_node_type{std::move(base_node)},
            _max_branch_size{max_branch_size},
            _remaining_branch_size{remaining_branch_size}
        {}

    public:

        explicit constexpr node_impl() = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            if (can_be_trimmed()) {
                return visit<true>(base_node_type::next_alt());
            } else {
                return std::nullopt;
            }
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            if (can_be_trimmed()) {
                return visit<false>(base_node_type::next_ref());
            } else {
                return std::nullopt;
            }
        }

        constexpr libjst::node_label_t<base_node_type> operator*() const noexcept {
            return *static_cast<base_node_type const &>(*this);
        }

    protected:

        constexpr breakpoint right_breakpoint() const {
            breakpoint original_bp = base_node_type::right_breakpoint();
            if (base_node_type::on_alternate_path()) {
                assert(static_cast<difference_type>(original_bp) + _remaining_branch_size > 0);
                return std::min<difference_type>(original_bp, static_cast<difference_type>(original_bp) + _remaining_branch_size);
            } else {
                return original_bp;
            }
        }
    private:

        constexpr bool can_be_trimmed() const noexcept {
            return _remaining_branch_size > 0;
        }

        template <bool is_alt_node, typename maybe_child_t>
        constexpr std::optional<node_impl> visit(maybe_child_t maybe_child) const {
            if (maybe_child) {
                if (is_alt_node && !base_node_type::on_alternate_path()) {
                    return branch_off_new(std::move(*maybe_child));
                } else if (base_node_type::on_alternate_path()) {
                    return branch_off_further<is_alt_node>(std::move(*maybe_child));
                } else { // nothing to spawn - remain in the reference branch.
                    return node_impl{std::move(*maybe_child), _max_branch_size, _remaining_branch_size};
                }
            } else {
                return std::nullopt;
            }
        }

        constexpr node_impl branch_off_new(base_node_type base_child) const noexcept {
            node_impl new_child{std::move(base_child), _max_branch_size, _remaining_branch_size};
            new_child._max_branch_size += breakpoint_span<true>();
            return new_child;
        }

        template <bool is_alt_node>
        constexpr node_impl branch_off_further(base_node_type base_child) const noexcept {
            node_impl new_child{std::move(base_child), _max_branch_size, _remaining_branch_size};
            new_child._remaining_branch_size -= new_child.breakpoint_span<is_alt_node>();
            return new_child;
        }

        template <bool is_alt_node>
        constexpr difference_type breakpoint_span() const noexcept {
            if constexpr (is_alt_node) {
                return std::ranges::ssize(libjst::alt_sequence(*base_node_type::left_variant()));
            } else {
                return base_node_type::right_breakpoint().value() - base_node_type::left_breakpoint().value();
            }
        }

        // constexpr bool on_variant_branch() const noexcept {
        //     return _on_variant_branch;
        // }

        // constexpr switch_to_variant_branch() const noexcept {
        //     _on_variant_branch = true;
        // }

        constexpr friend bool operator==(node_impl const & lhs, sink_impl const & rhs) noexcept {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    template <typename base_tree_t>
    class trim_tree_impl<base_tree_t>::sink_impl {
    private:
        friend trim_tree_impl;

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
        inline constexpr struct _trim
        {
            template <typename labelled_tree_t, std::unsigned_integral branch_size_t>
            constexpr auto operator()(labelled_tree_t && tree, branch_size_t const branch_size) const
                noexcept(std::is_nothrow_constructible_v<trim_tree_impl<std::remove_reference_t<labelled_tree_t>>>)
                -> trim_tree_impl<std::remove_reference_t<labelled_tree_t>>
            {
                using adapted_tree_t = trim_tree_impl<std::remove_reference_t<labelled_tree_t>>;
                return adapted_tree_t{(labelled_tree_t &&)tree, branch_size};
            }

            template <std::unsigned_integral branch_size_t>
            constexpr auto operator()(branch_size_t const branch_size) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, branch_size_t>)
                -> jst::contrib::closure_result_t<_trim, branch_size_t>
            {
                return jst::contrib::make_closure(_trim{}, branch_size);
            }
        } trim{};
    } // namespace _tree_adaptor

    using _tree_adaptor::trim;
}  // namespace libjst
