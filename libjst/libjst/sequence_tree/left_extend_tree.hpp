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
#include <libjst/sequence_tree/transform_tree.hpp>

namespace libjst
{

    template <typename base_label_t>
    class left_extended_sequence_label : public base_label_t {
    private:

        std::ptrdiff_t _left_extension{};
    public:

        using typename base_label_t::size_type;

        left_extended_sequence_label() = default;
        constexpr explicit left_extended_sequence_label(base_label_t base_label,
                                                        std::ptrdiff_t left_extension) noexcept :
            base_label_t{std::move(base_label)},
            _left_extension{left_extension}
        {}

        constexpr auto sequence(size_type left_pos = 0, size_type right_pos = base_label_t::npos) const {
            left_pos = std::max<std::ptrdiff_t>(0, left_pos - _left_extension);
            return base_label_t::sequence(left_pos, right_pos);
        }
    };

    class left_extend_factory {
    private:

        std::ptrdiff_t _left_extension{};
    public:

        left_extend_factory() = default;
        template <std::integral extension_t>
        constexpr explicit left_extend_factory(extension_t const left_extension) noexcept :
            _left_extension{static_cast<std::ptrdiff_t>(left_extension)}
        {}

        template <typename base_label_t>
        constexpr auto operator()(base_label_t && base_label) const noexcept {
            using wrapped_label_t = left_extended_sequence_label<std::remove_cvref_t<base_label_t>>;
            return wrapped_label_t{(base_label_t &&) base_label, _left_extension};
        }
    };

    namespace _tree_adaptor {
        inline constexpr struct _left_extend
        {
            template <typename labelled_tree_t, std::integral left_extension_t>
            constexpr auto operator()(labelled_tree_t && tree, left_extension_t left_extension) const
                -> decltype(transform((labelled_tree_t &&) tree, left_extend_factory{left_extension}))
                // noexcept(std::is_nothrow_constructible_v<
                //             transform_tree<std::remove_reference_t<labelled_tree_t>>, left_extension_t>)
                // -> transform_tree<std::remove_reference_t<labelled_tree_t>>
            {
                assert(left_extension >= 0);
                return transform((labelled_tree_t &&) tree, left_extend_factory{left_extension});
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
