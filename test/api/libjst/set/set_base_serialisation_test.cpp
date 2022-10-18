// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <sstream>

#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

#include <libcontrib/seqan/alphabet.hpp>

#include <libjst/set/concept_set.hpp>
#include <libjst/set/set_base.hpp>
#include <libjst/utility/bit_vector.hpp>
#include <libjst/variant/variant_snp.hpp>
#include <libjst/variant/variant_generic.hpp>
#include <libjst/variant/variant_store_composite.hpp>
#include <libjst/variant/variant_store_covered.hpp>

#include <libjst/set/serialiser_direct.hpp>
#include <libjst/set/serialiser_delegate.hpp>

namespace test {
template <typename other_t>
struct my_class
{
    other_t const & value;

private:
    // how do i get the wrappee stuff going on?
    template <typename archive_t>
    constexpr friend auto tag_invoke(std::tag_t<libjst::load>, my_class & me, archive_t & archive)
    {
        libjst::load_extern(archive, me.value); // load nested
        // archive(libjst::serialise_external{me.value});
    }

    template <typename archive_t>
    constexpr friend auto tag_invoke(std::tag_t<libjst::save>, my_class const & me, archive_t & archive)
    {
        // libjst::save_extern(archive, me.value)
        libjst::save_extern(archive, me.value);
    }
};

template <typename wrappee_t>
struct my_class_wrapper
{
    wrappee_t const & _wrappee; // now how do we get this loaded through?
    double _new_value{};

private:
    template <typename archive_t>
        requires (!std::same_as<std::remove_cvref_t<archive_t>, my_class_wrapper>)
    constexpr friend auto tag_invoke(std::tag_t<libjst::load>, my_class_wrapper & me, archive_t & archive)
    {
        libjst::load_extern(archive, me._wrappee); // call the nested object
        archive(me._new_value); // call its current object
    }

    template <typename archive_t>
        requires (!std::same_as<std::remove_cvref_t<archive_t>, my_class_wrapper>)
    constexpr friend auto tag_invoke(std::tag_t<libjst::save>, my_class_wrapper const & me, archive_t & archive)
    {
        // libjst::save_extern(archive, me._wrappee);
        libjst::save_extern(archive, me._wrappee);
        archive(me._new_value);
    }
};
}
// namespace libjst {



// TODO: later!
// namespace detail
// {
//     // TODO: what if the closure contains data?
//     template <typename jst_t, typename ...closures_t>
//         requires (sizeof...(closures_t) == 0)
//     struct compose_jst
//     {
//         using type = jst_t;
//     };

//     template <typename jst_t, typename next_closure_t, typename ...closures_t>
//     struct compose_jst<jst_t, next_closure_t, closures_t...>
//     {
//         using type = typename compose_jst<std::invoke_result_t<jst_t, next_closure_t>, ...closures_t>::type;
//     };

//     template <typename jst_model_t, typename ...closures_t>
//     using compose_jst_t = typename compose_jst<jst_model_t, closures_t...>::type;

// } // namespace detail

// template <typename jst_model_t, typename ...component_closures_t>
// class jst_serializer
// {
// private:

//     // we need to assemble the pipeline?
//     using jst_t = detail::compose_jst_t<jst_model_t, comonent_closures_t...>;

//     // what is he getting?
//     std::shared_ptr<jst_t const> _jst{}; // after handling it here we never remove const-ness?

// public:

//     template <typename >
//     jst_serializer()

//     // resetting the owened jst if any and sets the shared pointer to the external resource.
//     void reset(jst_t const &jst) noexcept
//     {
//         _jst.reset(std::addressof(jst));
//     }

//     jst_t const & get() const noexcept
//     {
//         return _jst.get(); // returns managed jst object.
//     }

//     // reference archive is provided from outside?
//     // maybe not?
//     template <typename archive_t>
//     void save(archive_t & archive)
//     {
//         assert(_jst != nullptr);
//         // now we need to get the assembly pipeline!
//         auto serialiser = assemble_serialiser_pipeline(archive,
//                                                        std::make_index_sequence<sizeof...(component_closures_t)>());
//         // now we assemble the pipeline

//     }

//     template <typename archive_t>
//     void load(archive_t & archive)
//     {

//     }

// private:

//     template <typename archive_t>
//     auto assemble_serialiser_pipeline(archive_t &archive) const
//     {
//         // we need the compnents and the memory location of them?
//         // so how do we get this done, if they are all wrapped?
//         auto assmble_impl = [&] <typename step_v, typename component_archive_t>
//         (step_v const &, component_archive_t && component_archive)
//         {
//             if constexpr (count_v::value == 0)
//                 return component_archive;
//             else
//                 return assemble_impl(std::integral_constant<size_t, step_v::value - 1>,
//                                      component_archive | delegate_serialiser()); // what element are we supposed to serialise here?
//         } (std::integral_constant<size_t, sizeof...(component_closures_t)>, archive) ;
//     }
// };
// } // namespace libjst

// how do we get these points automatically?
// get a serialiser -> does the assembly!
// get a director -> to select components to build and then create serialiser and build?

TEST(journaled_sequence_tree_serialiser_test, protoype)
{
    int data_out{10};
    int data_in{};

    test::my_class<int> ext_out{data_out};
    test::my_class_wrapper w1_out{ext_out, 0.5};
    test::my_class_wrapper w2_out{w1_out, 1.3};
    test::my_class_wrapper w3_out{w2_out, -24.1};

    test::my_class<int> ext_in{data_in};
    test::my_class_wrapper w1_in{ext_in, 0.0};
    test::my_class_wrapper w2_in{w1_in, 0.0};
    test::my_class_wrapper w3_in{w2_in, 0.0};

    // directly: build a archive pipeline

    // libjst::journaled_sequence_tree_serialiser serialiser_out{ext_out};
    // libjst::journaled_sequence_tree_serialiser serialiser_in{ext_in, my_loader};

    std::stringstream archive_stream{};
    {
        cereal::JSONOutputArchive output_archive(archive_stream);
        auto fst = libjst::direct_serialiser(data_out)
                 | libjst::delegate_serialiser(ext_out)
                 | libjst::delegate_serialiser(w1_out)
                 | libjst::delegate_serialiser(w2_out);
        // libjst::direct_archive arch{output_archive, data_out}; // now the archive for storing the external data.
        // libjst::delegate_archive arch2{arch, ext_out};
        // arch(ext_out); -> we could do a closure object here.
        auto arch = fst(output_archive);
        libjst::save(w3_out, arch);
    }

    std::cout << "stream: " << archive_stream.str() << "\n";
    EXPECT_NE(ext_in.value, ext_out.value);
    EXPECT_NE(data_in, data_out);
    EXPECT_NE(w1_in._new_value, w1_out._new_value);
    EXPECT_NE(w2_in._new_value, w2_out._new_value);
    EXPECT_NE(w3_in._new_value, w3_out._new_value);

    {
        cereal::JSONInputArchive input_archive(archive_stream);
        // libjst::direct_archive arch{input_archive, data_in}; // now the archive for storing the external data.
        auto fst = libjst::direct_serialiser(data_in)
                 | libjst::delegate_serialiser(ext_in)
                 | libjst::delegate_serialiser(w1_in)
                 | libjst::delegate_serialiser(w2_in);
        auto arch = fst(input_archive);
        libjst::load(w3_in, arch);
    }

    EXPECT_EQ(ext_in.value, ext_out.value);
    EXPECT_EQ(data_in, data_out);
    EXPECT_DOUBLE_EQ(w1_in._new_value, w1_out._new_value);
    EXPECT_DOUBLE_EQ(w2_in._new_value, w2_out._new_value);
    EXPECT_DOUBLE_EQ(w3_in._new_value, w3_out._new_value);
}

TEST(journaled_sequence_tree_serialiser_test, protoype_jst)
{
    using alphabet_t = jst::contrib::dna5;
    using sequence_t = std::vector<alphabet_t>;
    using snp_variant_t = libjst::snp_variant<alphabet_t>;
    using generic_variant_t = libjst::generic_variant<alphabet_t>;
    using coverage_t = libjst::bit_vector<>;

    using snp_store_t = std::vector<snp_variant_t>;
    using generic_store_t = std::vector<generic_variant_t>;
    using composite_store_t = libjst::variant_store_composite<snp_store_t, generic_store_t>;
    using covered_store_t = libjst::variant_store_covered<composite_store_t, libjst::bit_vector<>>;
    using value_t = std::ranges::range_value_t<covered_store_t>;

    using jst_t = libjst::set_base<sequence_t, covered_store_t>;

    std::vector<alphabet_t> base_sequence{seqan3::test::generate_sequence<alphabet_t>(200)};
    std::vector<alphabet_t> insertion_sequence{seqan3::test::generate_sequence<alphabet_t>(10)};

    snp_variant_t snp0{4, seqan3::assign_char_to('T', alphabet_t{})};
    snp_variant_t snp1{112, seqan3::assign_char_to('A', alphabet_t{})};
    generic_variant_t var0{44, insertion_sequence, 10};
    generic_variant_t var1{93, insertion_sequence, 0};
    generic_variant_t var2{154, {}, 1};

    jst_t jst_out{base_sequence, 4};

    EXPECT_TRUE((jst_out.insert(value_t{snp0, coverage_t{0, 0, 0, 1}})));
    EXPECT_TRUE((jst_out.insert(value_t{var0, coverage_t{0, 0, 1, 0}})));
    EXPECT_TRUE((jst_out.insert(value_t{var1, coverage_t{0, 1, 0, 0}})));
    EXPECT_TRUE((jst_out.insert(value_t{snp1, coverage_t{1, 0, 0, 0}})));
    EXPECT_TRUE((jst_out.insert(value_t{var2, coverage_t{0, 0, 1, 1}})));

    // test output stream and input stream in same buffer
    std::stringstream archive_stream{};
    {
        cereal::JSONOutputArchive output_archive(archive_stream);
        auto arch = output_archive | libjst::direct_serialiser(base_sequence);
        libjst::save(jst_out, arch);
    }

    std::vector<alphabet_t> base_sequence_in{};
    jst_t jst_in{base_sequence_in, 0}; // some default values.
    { // input arcihve
        cereal::JSONInputArchive input_archive(archive_stream);
        auto arch = input_archive | libjst::direct_serialiser(base_sequence_in);
        libjst::load(jst_in, arch);
    }

    EXPECT_RANGE_EQ(libjst::base_sequence(jst_in), libjst::base_sequence(jst_out));
    EXPECT_EQ(libjst::size(jst_in), libjst::size(jst_out));
    auto const & variant_store_out = libjst::variant_store(jst_out);
    auto const & variant_store_in = libjst::variant_store(jst_in);
    EXPECT_EQ(std::ranges::size(variant_store_in), std::ranges::size(variant_store_out));
    for (unsigned i = 0; i < std::ranges::size(variant_store_in); ++i)
    {
        EXPECT_EQ(libjst::position(variant_store_in[i]), libjst::position(variant_store_out[i]));
        EXPECT_EQ(libjst::deletion(variant_store_in[i]), libjst::deletion(variant_store_out[i]));
        EXPECT_RANGE_EQ(libjst::insertion(variant_store_in[i]), libjst::insertion(variant_store_out[i]));
        EXPECT_RANGE_EQ(libjst::coverage(variant_store_in[i]), libjst::coverage(variant_store_out[i]));
    }
}
