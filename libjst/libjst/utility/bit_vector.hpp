// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::utility::bit_vector.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/utility/container/aligned_allocator.hpp>
#include <libjst/utility/bit_vector_base.hpp>

namespace libjst
{

/*!\brief An allocator aware bit vector.
 *
 * \tparam allocator_t The type of the allocator to use.
 *
 * \details
 *
 * Implements an [allocator aware](https://en.cppreference.com/w/cpp/named_req/AllocatorAwareContainer) bit vector
 * on the basis of a std::vector with `uint64_t` as value type. The bit vector can be dynamically resized and
 * provides additional interfaces to apply efficient bit-operations on it.
 *
 * The reference type is a special proxy that provides access to a single bit. Note that it is not a real reference
 * but can be converter to a bool or assigned from a bool.
 */
template <typename allocator_t = std::allocator<uint64_t>>
//!\cond
    requires (!std::integral<allocator_t>)
//!\endcond
class bit_vector : public bit_vector_base<bit_vector<allocator_t>, allocator_t>
{
    //!\brief The crtp-base class.
    using base_t = bit_vector_base<bit_vector<allocator_t>, allocator_t>;

    //!\brief Befriend the base class to get access to the impl functions.
    friend base_t;

    //!\copydoc libjst::utility::bit_vector_base::chunk_type
    using typename base_t::chunk_type;

    static constexpr std::ptrdiff_t _unroll_factor{32}; //!< The number of chunks to roll out in a SIMD operation.

public:
    /*!\name Associated types
     * \{
     */
    //!\copydoc libjst::utility::bit_vector_base::iterator
    using typename base_t::iterator;
    //!\copydoc libjst::utility::bit_vector_base::const_iterator
    using typename base_t::const_iterator;
    //!\copydoc libjst::utility::bit_vector_base::value_type
    using typename base_t::value_type;
    //!\copydoc libjst::utility::bit_vector_base::reference
    using typename base_t::reference;
    //!\copydoc libjst::utility::bit_vector_base::const_reference
    using typename base_t::const_reference;
    //!\copydoc libjst::utility::bit_vector_base::size_type
    using typename base_t::size_type;
    //!\copydoc libjst::utility::bit_vector_base::difference_type
    using typename base_t::difference_type;
    //!\copydoc libjst::utility::bit_vector_base::allocator_type
    using typename base_t::allocator_type;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    bit_vector() = default; //!< Default.
    bit_vector(bit_vector const &) = default; //!< Default.
    bit_vector(bit_vector &&) = default; //!< Default.
    bit_vector & operator=(bit_vector const &) = default; //!< Default.
    bit_vector & operator=(bit_vector &&) = default; //!< Default.
    ~bit_vector() = default; //!< Default.

    /*!\brief Constructs the bit vector with `count` copies of elements with value `bit`.
     *
     * \param[in] count The number of elements to create the bit vector with.
     * \param[in] bit The bit to set during initialisation.
     * \param[in] alloc The allocator to use [optional].
     */
    constexpr bit_vector(size_type const count, bool const bit, allocator_t const & alloc = allocator_t{}) :
        base_t{alloc}
    {
        base_t::assign(count, bit);
    }

    /*!\brief Constructs the container initialised with the elements in `list`.
     *
     * \param[in] list An initialiser list with the bits set.
     * \param[in] alloc The allocator to use [optional].
     */
    constexpr bit_vector(std::initializer_list<bool> list, allocator_t const & alloc = allocator_t{})
        : base_t{alloc}
    {
        base_t::assign(list);
    }

    /*!\brief Constructs the container with `count` default-inserted instances of `bool`. No copies are made.
     *
     * \param[in] count The number of elements to create the bit vector with.
     * \param[in] alloc The allocator to use [optional].
     */
    constexpr bit_vector(size_type const count, allocator_t const & alloc = allocator_t{}) :
        bit_vector{count, bool{}, alloc}
    {}
    //!\}

    constexpr bool all() const noexcept
    {
        size_t data_size = base_t::as_base()->size(); // number of words
        for (size_t i = 0; i < data_size; ++i)
        {
            if (~base_t::as_base()->data()[i] != 0)
                return false;
        }
        return true;
    }

    constexpr bool any() const noexcept
    {
        size_t data_size = base_t::as_base()->size(); // number of words
        for (size_t i = 0; i < data_size; ++i)
        {
            if (base_t::as_base()->data()[i] > 0)
                return true;
        }
        return false;
    }

private:
    //!\brief Performs the binary bitwise-operation on the underlying chunks.
    template <typename binary_operator_t>
    static constexpr void binary_transform_impl(bit_vector & res,
                                                bit_vector const & lhs,
                                                bit_vector const & rhs,
                                                binary_operator_t && op) noexcept
    {
        assert(lhs.size() == rhs.size());
        assert(res.size() == lhs.size());

        size_t unroll_offset = _unroll_factor;
        size_t data_size = lhs.as_base()->size(); // number of words

        for (; unroll_offset < data_size; unroll_offset += _unroll_factor)
        {
            auto start = unroll_offset - _unroll_factor;
            for (size_t i = 0; i < _unroll_factor; ++i)
            {
                res.as_base()->data()[start + i] = op(lhs.as_base()->data()[start + i], rhs.as_base()->data()[start + i]);
            }
        }

        for (size_t i = unroll_offset - _unroll_factor; i < data_size; ++i)
        {
            res.as_base()->data()[i] = op(lhs.as_base()->data()[i], rhs.as_base()->data()[i]);
        }
    }

    //!\brief Performs the binary bitwise-operation on the underlying chunks.
    template <typename binary_operator_t>
    static constexpr void unary_transform_impl(bit_vector & res,
                                               bit_vector const & lhs,
                                               binary_operator_t && op) noexcept
    {
        assert(res.size() == lhs.size());

        size_t unroll_offset = _unroll_factor;
        size_t data_size = lhs.as_base()->size(); // number of words

        for (; unroll_offset < data_size; unroll_offset += _unroll_factor)
        {
            auto start = unroll_offset - _unroll_factor;
            for (size_t i = 0; i < _unroll_factor; ++i)
            {
                res.as_base()->data()[start + i] = op(lhs.as_base()->data()[start + i]);
            }
        }

        for (size_t i = unroll_offset - _unroll_factor; i < data_size; ++i)
        {
            res.as_base()->data()[i] = op(lhs.as_base()->data()[i]);
        }
    }

    //!\brief Computes the minimal size needed for the host vector.
    //!\param[in] count The number of bits to allocate memory for.
    constexpr size_type host_size_impl(size_type const count) const noexcept
    {
        return static_cast<uint64_t>(base_t::chunks_needed(count));
    }
};
}  // namespace libjst::utility

//!\cond
namespace cereal
{
template <typename archive_t, typename allocator_t>
struct specialize<archive_t, libjst::bit_vector<allocator_t>, cereal::specialization::member_load_save>
{};
} // namespace cereal
//!\endcond
