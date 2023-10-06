// Adds basic serialisation concepts and default functions using tag_invoke based CPOs

#pragma once

#include <libjst/utility/tag_invoke.hpp>

#if __has_include(<cereal/cereal.hpp>)
#include <cereal/cereal.hpp>

namespace libjst::detail {

    inline constexpr bool has_cereal = true;

    template <typename archive_t>
    concept cereal_input_archive = std::derived_from<archive_t, cereal::detail::InputArchiveBase>;

    template <typename archive_t>
    concept cereal_output_archive = std::derived_from<archive_t, cereal::detail::OutputArchiveBase>;

}

#else

namespace libjst::detail {

    inline constexpr bool has_cereal = false;

    template <typename archive_t>
    concept cereal_input_archive = true;

    template <typename archive_t>
    concept cereal_output_archive = true;
}

#endif

/**
 * @defgroup serialisation Serialisation
 *
 * General serialisation support.
 */

namespace libjst {

    namespace _load {

        inline constexpr struct _cpo {
            /**
             * @brief Calling this CPO will load an object from a given input archive.
             *
             * @tparam object_t The object type.
             * @tparam iarchive_t The input archive type.
             *                    Must be an <a href="https://uscilab.github.io/cereal/assets/doxygen/classcereal_1_1InputArchive.html">cereal::InputArchive</a>.
             * @param object The object to load.
             * @param iarchive The input archive to load from.
             *
             * To use this CPO, the object type must offer one of the following three function signatures:
             *   - A member function with the signature `void load(iarchive_t&)`;
             *   - A free function with the signature `void load(object_t&, iarchive_t&)`;
             *   - Or a libjst::tag_invocable overload with the signature `friend void tag_invoke(libjst::tag_t<libjst::load>, object_t&, iarchive_t&)`.
             *
             * @ingroup serialisation
             */
            template <typename object_t, detail::cereal_input_archive iarchive_t>
                requires (detail::has_cereal) && (libjst::tag_invocable<_cpo, object_t&, iarchive_t&>)
            constexpr auto operator()(object_t& object, iarchive_t& iarchive) const
                noexcept(libjst::is_nothrow_tag_invocable_v<_cpo, object_t&, iarchive_t &>)
                -> libjst::tag_invoke_result_t<_cpo, object_t&, iarchive_t&> {
                return libjst::tag_invoke(_cpo{}, object, iarchive);
            }

            /**
             * @brief Default overload if cereal is not available.
             *
             * This overload will always throw a std::runtime_error.
             */
            template <typename object_t, typename iarchive_t>
                requires (!detail::has_cereal)
            constexpr void operator()(object_t&, iarchive_t&) const {
                throw std::runtime_error("libjst::load: cereal is not available");
            }

        private:
            /**
             * @brief This overload is only enabled if the member load function is available.
             */
            template <typename object_t, typename iarchive_t>
                requires (requires(object_t& object, iarchive_t& iarchive) {
                            { object.load(iarchive) } -> std::same_as<void>;
                         })
            friend void tag_invoke(_cpo, object_t& object, iarchive_t& iarchive)
                noexcept(noexcept(object.load(iarchive))) {
                object.load(iarchive);
            }

            /**
             * @brief This overload is only enabled if the free load function is available
             *        and the member load function is not available
             */
            template <typename object_t, typename iarchive_t>
                requires (!requires(object_t& object, iarchive_t& iarchive) {
                            { object.load(iarchive) } -> std::same_as<void>;
                         }) &&
                         (requires(object_t& object, iarchive_t& iarchive) {
                            { load(object, iarchive) } -> std::same_as<void>;
                         })
            friend void tag_invoke(_cpo, object_t& object, iarchive_t& iarchive)
                noexcept(noexcept(load(object, iarchive))) {
                load(object, iarchive);
            }
        } load;
    } // namespace _load

    using _load::load;

    // Add CPO implementation to save objects using cereal output archives
    namespace _save
    {
        inline constexpr struct _cpo {
            /**
             * @brief Calling this CPO will save an object to a given output archive.
             *
             * @tparam object_t The object type.
             * @tparam oarchive_t The output archive type.
             *                    Must be an <a href="https://uscilab.github.io/cereal/assets/doxygen/classcereal_1_1OutputArchive.html">cereal::OutputArchive</a>.
             *
             * @param object The object to save.
             * @param oarchive The output archive to save to.
             *
             * To use this CPO, the object type must offer one of the following three function signatures:
             *  - A member function with the signature `void save(oarchive_t&) const`;
             *  - A free function with the signature `void save(object_t const &, oarchive_t&)`;
             *  - Or a libjst::tag_invocable overload with the signature `friend void tag_invoke(libjst::tag_t<libjst::save>, object_t const &, oarchive_t&)`.
             *
             * @ingroup serialisation
             */
            template <typename object_t, detail::cereal_output_archive oarchive_t>
                requires (detail::has_cereal) && (libjst::tag_invocable<_cpo, object_t&, oarchive_t&>)
            constexpr auto operator()(object_t const & object, oarchive_t& oarchive) const
                noexcept(libjst::is_nothrow_tag_invocable_v<_cpo, object_t const &, oarchive_t &>)
                -> libjst::tag_invoke_result_t<_cpo, object_t const &, oarchive_t&> {
                return libjst::tag_invoke(_cpo{}, object, oarchive);
            }

            /**
             * @brief This overload is only enabled if cereal is not available.
             *
             * This overload will always throw a std::runtime_error.
             */
            template <typename object_t, typename oarchive_t>
                requires (!detail::has_cereal)
            void operator()(object_t&, oarchive_t&) const {
                throw std::runtime_error("libjst::save: cereal is not available");
            }

        private:

            /**
             * @brief Default overload using the save member function.
             */
            template <typename object_t, typename oarchive_t>
                requires (requires(object_t const & object, oarchive_t& oarchive) {
                            { object.save(oarchive) } -> std::same_as<void>;
                         })
            friend void tag_invoke(_cpo, object_t const & object, oarchive_t& oarchive)
                noexcept(noexcept(object.save(oarchive))) {
                object.save(oarchive);
            }

            /**
             * @brief Default overload using a free save function.
             */
            template <typename object_t, typename oarchive_t>
                requires (!requires(object_t const & object, oarchive_t& oarchive) {
                            { object.save(oarchive) } -> std::same_as<void>;
                         }) &&
                         (requires(object_t const & object, oarchive_t& oarchive) {
                            { save(object, oarchive) } -> std::same_as<void>;
                         })
            friend void tag_invoke(_cpo, object_t const & object, oarchive_t& oarchive)
                noexcept(noexcept(save(object, oarchive))) {
                save(object, oarchive);
            }

        } save;
    }

    using _save::save;
}
