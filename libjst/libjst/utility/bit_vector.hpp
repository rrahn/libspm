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
template <typename allocator_t = std::allocator<bool>>
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

private:
    //!\brief Performs the binary bitwise-operation on the underlying chunks.
    template <typename binary_operator_t>
    constexpr bit_vector & binary_transform_impl(bit_vector const & rhs, binary_operator_t && op) noexcept
    {
        std::ranges::transform(*base_t::as_base(), *rhs.as_base(), base_t::as_base()->data(),
                               [&] (chunk_type const & left_chunk, chunk_type const & right_chunk) -> chunk_type
        {
            return op(left_chunk, right_chunk);
        });

        return *this;
    }

    //!\brief Ensures that the base bit vector has enough memory.
    constexpr void reserve_impl(size_type const count) noexcept
    {
        base_t::as_base()->reserve(base_t::chunks_needed(count));
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
