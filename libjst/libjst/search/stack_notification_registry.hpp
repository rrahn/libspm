// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::search_stack_notification_registry.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <memory>
#include <vector>

#include <libjst/search/stack_observer_concept.hpp>

namespace libjst
{

/*!\brief A stack notification registry for observable search algorithms.
 *
 * \details
 *
 * This registry manages a set of libjst::search_stack_observer's. During construction a set of observers can register
 * in the registry. Whenever the registry is notified about a stack event (push/pop), it will notify the registered
 * observers. It uses some internal runtime polymorphism which is not exposed observers. Each observer must model
 * the libjst::search_stack_observer concept.
 * The attached observers must outlive the registry, otherwise calling a notification event on a dangling reference
 * will lead to undefined behaviour.
 */
class search_stack_notification_registry
{
private:

    //!\brief The abstract base class to implement the observable concept.
    struct observable
    {
        observable() = default; //!< Default.
        virtual ~observable() = default; //!< Virtual and default.

        //!\brief Notifies the attached observer about a pop event.
        virtual void notify_pop() = 0;
        //!\brief Notifies the attached observer about a push event.
        virtual void notify_push() = 0;
    };

    /*!\brief A custom observer specialised over the user provided stack observer.
     *
     * \tparam observer_t The observer type to be stored.
     *
     * \details
     *
     * Handles the attached observer and forwards the notification events to the attached observers.
     */
    template <search_stack_observer observer_t>
    struct observer_handle : public observable
    {
        observer_handle() = default; //!< Default.

        //!\brief Constructs a new observer handle to manage the attached observer.
        //!\param[in] observer A reference to the attached observer.
        observer_handle(observer_t & observer) : _handle{std::addressof(observer)}
        {}

        //!\copydoc libjst::search_stack_notification_registry::observable::notify_pop()
        void notify_pop() override
        {
            assert(_handle != nullptr);
            _handle->on_pop();
        }

        //!\copydoc libjst::search_stack_notification_registry::observable::notify_push()
        void notify_push() override
        {
            assert(_handle != nullptr);
            _handle->on_push();
        }

        observer_t * _handle{nullptr}; //!\< The handle to the attached observer.
    };

    std::vector<std::unique_ptr<observable>> _observer_list; //!< The list of attached observers.

public:

    /*!\name Constructors, destructor, and assignment
     * \{
     */
    search_stack_notification_registry() = default; //!< Default.
    search_stack_notification_registry(search_stack_notification_registry const &) = delete; //!< Deleted.
    search_stack_notification_registry(search_stack_notification_registry &) = default; //!< Default.
    search_stack_notification_registry & operator=(search_stack_notification_registry const &) = delete; //!< Deleted.
    search_stack_notification_registry & operator=(search_stack_notification_registry &) = default; //!< Default.
    ~search_stack_notification_registry() = default; //!< Default.

    //!\brief Constructs a new registry with a set of observers to attach.
    //!\tparam observer_t Template parameter pack for the observer types; must model libjst::search_stack_observer.
    template <search_stack_observer ...observer_t>
    //!\cond
        requires (sizeof...(observer_t) >= 1)
    //!\endcond
    explicit search_stack_notification_registry(observer_t & ...observer)
    {
        (_observer_list.emplace_back(std::make_unique<observer_handle<observer_t>>(observer)), ...);
    }
    //!\}

protected:
    //!\brief Notifies the attached observer about a push event.
    constexpr void notify_push()
    {
        notify_all([] (std::unique_ptr<observable> & observer_ptr) { observer_ptr->notify_push(); } );
    }

    //!\brief Notifies the attached observer about a pop event.
    constexpr void notify_pop()
    {
        notify_all([] (std::unique_ptr<observable> & observer_ptr) { observer_ptr->notify_pop(); } );
    }

private:
    //!\brief Notifies all attached observers about the incoming stack event.
    template <typename notification_event_t>
    constexpr void notify_all(notification_event_t && event)
    {
        std::ranges::for_each(_observer_list, event);
    }
};

}  // namespace libjst
