// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::journaled_sequence_tree partitioned.
 * \author Tom Lukas Lankenau <tom.lankenau AT fu-berlin.de>
 */

#pragma once

#include <vector>

#include <seqan3/core/detail/strong_type.hpp>

#include <libjst/journal_sequence_tree_context_enumerator.hpp>
#include <libjst/journal_sequence_tree_range_agent.hpp>
#include <libjst/journaled_sequence_tree.hpp>

namespace libjst
{

//!\brief A strong type to pass a context size object.
struct context_size : public seqan3::detail::strong_type<uint32_t, context_size>
{
    //!\brief The base type.
    using base_t = seqan3::detail::strong_type<uint32_t, context_size>;
    using base_t::base_t;
};

//!\brief A strong type to pass a bin index object.
struct bin_index : public seqan3::detail::strong_type<uint32_t, bin_index>
{
    //!\brief The base type.
    using base_t = seqan3::detail::strong_type<uint32_t, bin_index>;
    using base_t::base_t;
};

/*!\brief A partitioned journal sequence tree.
 *
 * \tparam jst_t The type of the journal sequence tree to wrap.
 *
 * \details
 *
 * This wrapper handles a collection of libjst::detail::journal_sequence_tree_traverser_model with non-overlapping
 * intervals over the given jst. Later the jst can be traversed in bins by constructing the respective agent for a
 * particular bin. The partitioned jst is a wrapper class and cannot be default initialised.
 */
template <typename jst_t>
class journal_sequence_tree_partitioned
{
public:
    //!\brief The type of the traverser model to manage.
    using traverser_model_t = detail::journal_sequence_tree_traverser_model<jst_t>;
    //!\brief The type of the context enumerator to use.
    using context_enumerator_type = typename jst_t::context_enumerator_type;
    //!\brief The type of the range agent to use.
    using range_agent_type = typename jst_t::range_agent_type;
    //!\brief The position type.
    using position_type = typename jst_t::position_type;

private:
    std::vector<traverser_model_t> _bins{}; //!\< The container stroing the model for each bin.
    jst_t const * _jst{}; //!\< A shared pointer to the underlying jst.

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr journal_sequence_tree_partitioned() = delete; //!< Delete.
    constexpr journal_sequence_tree_partitioned(journal_sequence_tree_partitioned const &) = default; //!< Default.
    constexpr journal_sequence_tree_partitioned(journal_sequence_tree_partitioned &&) = default; //!< Default.
    constexpr journal_sequence_tree_partitioned & operator=(journal_sequence_tree_partitioned const &)
        = default; //!< Default.
    constexpr journal_sequence_tree_partitioned & operator=(journal_sequence_tree_partitioned &&)
        = default; //!< Default.
    ~journal_sequence_tree_partitioned() = default; //!< Default.

    /*!\brief Constructs the partitioned jst from a jst pointer and a bin count.
     *
     * \param[in] jst A pointer to the jst to wrap.
     * \param[in] bin_size The size of a single bin. Defaults to infinity and must be `> 0`.
     */
    explicit journal_sequence_tree_partitioned(jst_t const * jst,
                                               size_t bin_size = std::numeric_limits<size_t>::max()) : _jst{jst}
    {
       assert(_jst != nullptr);
       assert(bin_size > 0);

        for (size_t ref_idx = 0; ref_idx < _jst->reference().size(); ++ref_idx)
        {
            // Invariant: bin_count >= 1
            size_t bin_count = std::max<size_t>(1, (_jst->reference_at(ref_idx).size() + bin_size - 1) / bin_size);

            for (size_t i = 0; i < bin_count; i++)
                _bins.push_back(traverser_model_t{_jst,
                                                  position_type{ref_idx, i * bin_size},
                                                  position_type{ref_idx, (i + 1) * bin_size + bin_overlap}});
        }
    }
    //!\}

    //!\brief Returns a copy of the bin at the given bin index.
    traverser_model_t bin_at(size_t bin_index) const noexcept
    {
        assert(bin_index < _bins.size());
        return _bins[bin_index];
    }

    //!\brief Returns the bin count.
    size_t bin_count() const noexcept
    {
        return _bins.size();
    }

    //!\brief Returns a new context enumerator from the given bin and the context size.
    context_enumerator_type context_enumerator(libjst::context_size const context_size,
                                               libjst::bin_index const bin_index) const
    {
        check_valid_bin_index(bin_index);
        return context_enumerator_type{_bins[bin_index.get()], context_size.get()};
    }

    //!\brief Returns a new context enumerator from the given bin and the context size.
    template <search_stack_observer ...observer_t>
    range_agent_type range_agent(libjst::context_size const context_size,
                                 libjst::bin_index const bin_index,
                                  observer_t & ...observer) const
    {
        check_valid_bin_index(bin_index);
        return range_agent_type{_bins[bin_index.get()], context_size.get(), observer...};
    }

    //!\brief Returns the sequence positions at the given coordinate.
    auto sequence_positions_at(journal_sequence_tree_coordinate const & coordinate) const
        -> decltype(_jst->sequence_positions_at(coordinate))
    {
        return _jst->sequence_positions_at(coordinate);
    }

    /*!\name Serialisation
     * \{
     */
    /*!\brief Saves this partitioned jst to the given output archive.
     *
     * \tparam output_archive_t The type of the output_archive; must model seqan3::cereal_output_archive.
     *
     * \param[in, out] archive The archive to serialise this object to.
     *
     * \details
     *
     * Manually saves the serialised traverser models to the archive. First the number of contained models is
     * stored. Subsequently, all models are stored consecutively.
     */
    template <seqan3::cereal_output_archive output_archive_t>
    void save(output_archive_t & archive) const
    {
        assert(_jst != nullptr);
        archive(_bins.size());
        std::ranges::for_each(_bins, [&] (traverser_model_t const & model)
        {
            model.save(archive);
        });
    }

    /*!\brief Loads this partitioned jst from the given input archive.
     *
     * \tparam input_archive_t The type of the input_archive; must model seqan3::cereal_input_archive.
     *
     * \param[in, out] archive The archive to serialise this object from.
     *
     *
     * details
     *
     * Manually loads the serialised traverser models from the archive. First the number of contained models is
     * loaded. Subsequently, all models are loaded in order from the archive initialised with the associated jst.
     */
    template <seqan3::cereal_input_archive input_archive_t>
    void load(input_archive_t & archive)
    {
        assert(_jst != nullptr);
        size_t bin_size{};
        archive(bin_size);
        _bins.resize(bin_size);
        std::ranges::for_each(_bins, [&] (traverser_model_t & model)
        {
            model.load(archive, _jst);
        });
    }
    //!\}

private:
    //!\brief Checks if the given bin index is valid.
    void check_valid_bin_index(libjst::bin_index const bin_index) const
    {
        if (bin_index.get() >= bin_count())
            throw std::out_of_range{"The bin index: " + std::to_string(bin_index.get()) + " is out of range!"};
    }
};

} // namespace libjst
