// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::utility::bit_vector_simd.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/utility/simd/algorithm.hpp>
#include <seqan3/utility/simd/simd.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>

#include <libjst/utility/bit_vector_base.hpp>

namespace libjst
{

/*!\brief A vectorised version of the bit vector.
 *
 * \tparam allocator_t The type of the allocator to use. Defaults to seqan3::aligned_allocator.
 *
 * \details
 *
 * Implements a bit vector with vectorised operations binary operations. The host memory should be aligned to a
 * multiple of the simd vector alignment specification. It always reserves more memory in order to read a full
 * vector from the underlying host memory without accessing memory that was not reserved by the bit vector.
 *
 * The reference type is a special proxy that provides access to a single bit. Note that it is not a real reference
 * but can be converter to a bool or assigned from a bool.
 */
template <typename allocator_t = seqan3::aligned_allocator<bool, alignof(seqan3::simd_type_t<uint64_t>)>>
//!\cond
    requires (!std::integral<allocator_t>)
//!\endcond
class bit_vector_simd : public bit_vector_base<bit_vector_simd<allocator_t>, allocator_t>
{
    //!\brief The crtp-base class.
    using base_t = bit_vector_base<bit_vector_simd<allocator_t>, allocator_t>;

    //!\brief Befriend the base class to get access to the impl functions.
    friend base_t;

    //!\copydoc libjst::utility::bit_vector_base::chunk_type
    using typename base_t::chunk_type;
    //!\brief The simd type for the chunk.
    using simd_type = seqan3::simd_type_t<chunk_type>;

    //!\brief How many chunks can be stored in one vector.
    static constexpr size_t chunks_per_vector = seqan3::simd_traits<simd_type>::length;

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
    bit_vector_simd() = default; //!< Default.
    //!\brief Copy constructor.
    bit_vector_simd(bit_vector_simd const & other) : base_t{}
    {
        *this = other;
    }

    bit_vector_simd(bit_vector_simd &&) = default; //!< Default.
    //!\brief Copy assignment.
    bit_vector_simd & operator=(bit_vector_simd const & other)
    {
        base_t::assign(other.begin(), other.end());
        return *this;
    }

    bit_vector_simd & operator=(bit_vector_simd &&) = default; //!< Default.
    ~bit_vector_simd() = default; //!< Default.

    /*!\brief Constructs the bit vector with `count` copies of elements with value `bit`.
     *
     * \param[in] count The number of elements to create the bit vector with.
     * \param[in] bit The bit to set during initialisation.
     * \param[in] alloc The allocator to use [optional].
     */
    constexpr bit_vector_simd(size_type const count, bool const bit, allocator_t const & alloc = allocator_t{}) :
        base_t{alloc}
    {
        base_t::assign(count, bit);
    }

    /*!\brief Constructs the container initialised with the elements in `list`.
     *
     * \param[in] list An initialiser list with the bits set.
     * \param[in] alloc The allocator to use [optional].
     */
    constexpr bit_vector_simd(std::initializer_list<bool> list, allocator_t const & alloc = allocator_t{})
        : base_t{alloc}
    {
        base_t::assign(list);
    }

    /*!\brief Constructs the container with `count` default-inserted instances of `bool`. No copies are made.
     *
     * \param[in] count The number of elements to create the bit vector with.
     * \param[in] alloc The allocator to use [optional].
     */
    constexpr bit_vector_simd(size_type const count, allocator_t const & alloc = allocator_t{}) :
        bit_vector_simd{count, bool{}, alloc}
    {}
    //!\}

    //!\brief Performs binary NOT.
    constexpr bit_vector_simd operator~() const noexcept
    {
        bit_vector_simd tmp(base_t::size());

        chunk_type const * lhs_data = base_t::as_base()->data();
        chunk_type const * end = base_t::as_base()->data() + min_capacity(base_t::size());
        chunk_type * rhs_data = tmp.as_base()->data();

        for (; lhs_data != end; lhs_data += chunks_per_vector, rhs_data += chunks_per_vector)
            seqan3::store(rhs_data, ~seqan3::load<simd_type>(lhs_data));

        return tmp;
    }

private:
    //!\brief Performs the binary bitwise-operation on the underlying chunks.
    template <typename binary_operator_t>
    constexpr bit_vector_simd & binary_transform_impl(bit_vector_simd const & rhs, binary_operator_t && op) noexcept
    {
        chunk_type * lhs_data = base_t::as_base()->data();
        chunk_type const * end = base_t::as_base()->data() + min_capacity(base_t::size());
        chunk_type const * rhs_data = rhs.as_base()->data();

        for (; lhs_data != end; lhs_data += chunks_per_vector, rhs_data += chunks_per_vector)
            seqan3::store(lhs_data, op(seqan3::load<simd_type>(lhs_data), seqan3::load<simd_type>(rhs_data)));

        return *this;
    }

    //!\brief Ensures that the base bit vector has enough memory.
    constexpr void reserve_impl(size_type const count) noexcept
    {
        base_t::as_base()->reserve(min_capacity(count));
    }

    /*!\brief Computes the minimal capacity needed to be allocated for the underlying vector.
     *
     * \param[in] count The number of bits to reserve memory for.
     *
     * \details
     *
     * In the simd approach we guarantee that the underlying vector allocates enough memory to always read an entire
     * simd vector without reading from unallocated memory.
     * It is not important that the memory is properly initialised as the overlapping data is ignored in any case.
     */
    constexpr size_t min_capacity(size_type const count) const noexcept
    {
        return (base_t::chunks_needed(count) + chunks_per_vector - 1) & -chunks_per_vector;
    }
};

}  // namespace libjst::utility
