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
#include <libcontrib/copyable_box.hpp>

#include <libjst/sequence_tree/concept.hpp>
#include <libjst/sequence_tree/journaled_sequence_label.hpp>
#include <libjst/variant/concept.hpp>

namespace libjst
{

    template <typename wrapped_tree_t>
    class labelled_tree {
    private:
        using base_node_type = libjst::tree_node_t<wrapped_tree_t>;
        using sink_type = libjst::tree_sink_t<wrapped_tree_t>;
        using base_cargo_type = libjst::tree_label_t<wrapped_tree_t>;

        using boundary_type = typename base_node_type::position_type;
        using delta_reference = typename boundary_type::delta_reference;
        using position_type = libjst::variant_position_t<delta_reference>;
        using sequence_type = libjst::alt_sequence_t<delta_reference>;

        // static_assert(std::same_as<void, sequence_type>);

        using label_strategy_type = journaled_sequence_label<position_type, sequence_type>;
        using tree_box_t = jst::contrib::copyable_box<wrapped_tree_t>;

        class node_impl;
        class cargo_impl;

        tree_box_t _wrappee{};

    public:

        template <typename wrappee_t>
            requires (!std::same_as<std::remove_cvref_t<wrappee_t>, labelled_tree> &&
                      std::constructible_from<wrapped_tree_t, wrappee_t>)
        constexpr explicit labelled_tree(wrappee_t && wrappee) noexcept : _wrappee{(wrappee_t &&)wrappee}
        {
        }

        constexpr node_impl root() const noexcept {
            label_strategy_type initial_label{data().source()};
            base_node_type root_base = libjst::root(_wrappee.value());
            initial_label.reset_positions(libjst::position(root_base.low_boundary()),
                                          libjst::position(root_base.high_boundary()));
            return node_impl{std::move(root_base), std::move(initial_label)};
        }

        constexpr sink_type sink() const noexcept {
            return libjst::sink(_wrappee.value());
        }

        constexpr auto const & data() const noexcept {
            return _wrappee->data();
        }
    };

    template <typename wrapped_tree_t>
    class labelled_tree<wrapped_tree_t>::node_impl : public base_node_type {
    private:

        friend labelled_tree;

        using base_t = base_node_type;

        label_strategy_type _label{};

        explicit constexpr node_impl(base_t && base_node, label_strategy_type label) noexcept :
                base_t{std::move(base_node)},
                _label{std::move(label)}
        {}

    public:

        node_impl() = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            return visit_next<true>(base_t::next_alt());
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            return visit_next<false>(base_t::next_ref());
        }

        constexpr cargo_impl operator*() const noexcept {
            return cargo_impl{this};
        }

    private:

        template <bool is_alt>
        constexpr std::optional<node_impl> visit_next(std::optional<base_t> && base_child) const noexcept {
            if (base_child.has_value()) {
                // I need to check between alt_sequence and sequence between interval!
                label_strategy_type child_label{_label};
                if constexpr (is_alt) {
                    child_label.record_variant(*(base_child->low_boundary()));
                }
                return node_impl{std::move(*base_child), std::move(child_label)};
            }
            return std::nullopt;
        }

        constexpr friend bool operator==(node_impl const & lhs, sink_type const & rhs) noexcept
        {
            return static_cast<base_t const &>(lhs) == rhs;
        }
    };

    template <typename wrapped_tree_t>
    class labelled_tree<wrapped_tree_t>::cargo_impl : public base_cargo_type {
    private:
        [[no_unique_address]] node_impl const * _node{};

        friend labelled_tree;

        constexpr explicit cargo_impl(node_impl const * node) :
            base_cargo_type{*(static_cast<base_node_type const &>(*node))},
            _node{node}
        {}
    public:

        using size_type = typename label_strategy_type::size_type;

        static constexpr size_type npos{label_strategy_type::npos};

        cargo_impl() = default;

        constexpr auto sequence() const noexcept {
            return sequence(libjst::position(_node->low_boundary()), libjst::position(_node->high_boundary()));
        }

    protected:

        constexpr auto sequence(size_type const first, size_type const last) const noexcept {
            assert(_node != nullptr);
            assert(first <= last);
            return _node->_label.sequence(first, last);
        }
    };


    // template <typename derived_t, typename extended_node_t, typename label_kind_t>
    // class node_label_extension {
    // private:

    //     friend derived_t;

    //     using breakend_iterator = std::remove_cvref_t<decltype(std::declval<extended_node_t const &>().low_boundary().get_breakend())>;
    //     using variant_type = std::iter_reference_t<breakend_iterator>;
    //     using sequence_type = libjst::alt_sequence_t<variant_type>;
    //     using position_type = libjst::variant_position_t<variant_type>;

    //     using label_strategy_type = journaled_sequence_label<position_type, sequence_type>;

    //     template <typename base_label_t>
    //     class label_impl;

    //     label_strategy_type _label_strategy{}; // the label information

    // private:

    //     node_label_extension() = default;

    //     constexpr void initialise() {
    //         _label_strategy = label_strategy_type(as_derived(*this).rcs_store().source());
    //         _label_strategy.reset_positions(as_derived(*this).low_breakend(),
    //                                         as_derived(*this).high_breakend());
    //     };

    //     constexpr node_label_extension notify(extended_node_t const & child_node) const {
    //         node_label_extension child_extension{*this}; // copy current extension into child.
    //         if (as_derived(child_node).is_alt_node()) {
    //             child_extension._label_strategy.record_variant(as_derived(child_node).left_variant()); //TODO
    //         } else {
    //             // now we could be in a multibranch state!
    //             // thus we need to set the correct left position when we go into state A!

    //             auto right_bp = as_derived(child_node).high_breakend();
    //             auto left_bp = as_derived(child_node).low_breakend();
    //             child_extension._label_strategy.reset_positions(std::min(left_bp, right_bp), right_bp);
    //         }
    //         return child_extension;
    //     }

    //     static constexpr derived_t & as_derived(node_label_extension & me) noexcept {
    //         return static_cast<derived_t &>(me);
    //     }

    //     static constexpr derived_t & as_derived(extended_node_t & me) noexcept {
    //         return static_cast<derived_t &>(me);
    //     }

    //     static constexpr derived_t const & as_derived(node_label_extension const & me) noexcept {
    //         return static_cast<derived_t const &>(me);
    //     }

    //     static constexpr derived_t const & as_derived(extended_node_t const & me) noexcept {
    //         return static_cast<derived_t const &>(me);
    //     }

    // public:

    //     template <typename base_label_t>
    //     constexpr auto label(base_label_t && base_label) const noexcept {
    //         return label_impl<std::remove_cvref_t<base_label_t>>{(base_label_t &&) base_label,
    //                                                              std::addressof(_label_strategy)};
    //     }
    // };

    // template <typename derived_t, typename extended_node_t, typename label_kind_t>
    // template <typename base_label_t>
    // class node_label_extension<derived_t, extended_node_t, label_kind_t>::label_impl : public base_label_t {
    // private:

    //     [[no_unique_address]] label_strategy_type const * _label_strategy{};

    //     friend node_label_extension;

    //     constexpr explicit label_impl(base_label_t base_label, label_strategy_type const * strategy) :
    //         base_label_t{std::move(base_label)},
    //         _label_strategy{strategy}
    //     {}
    // public:

    //     using size_type = typename label_strategy_type::size_type;

    //     static constexpr size_type npos{label_strategy_type::npos};

    //     label_impl() = default;

    //     constexpr auto sequence(size_type const first = 0, size_type const last = label_strategy_type::npos) const noexcept {
    //         assert(_label_strategy != nullptr);
    //         assert(first <= last);
    //         return _label_strategy->sequence(first, last);
    //     }
    // };

    namespace _tree_adaptor {
        struct _labelled {

            template <typename tree_t, typename ...args_t>
            constexpr auto operator()(tree_t && tree, args_t &&... args) const
                noexcept(std::is_nothrow_constructible_v<labelled_tree<std::remove_reference_t<tree_t>>>)
                -> labelled_tree<std::remove_reference_t<tree_t>>
            {
                using adapted_tree_t = labelled_tree<std::remove_reference_t<tree_t>>;
                return adapted_tree_t{(tree_t &&)tree, (args_t &&)args...};
            }

            template <typename ...args_t>
            constexpr auto operator()(args_t &&... args) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, _labelled,  args_t...>)
                -> jst::contrib::closure_result_t<_labelled, args_t...>
            {
                return jst::contrib::make_closure(_labelled{}, (args_t &&)args...);
            }
        };

        inline constexpr _labelled labelled{};
    } // namespace _tree_adaptor

    using _tree_adaptor::labelled;
}  // namespace libjst
