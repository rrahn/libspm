// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a stack observer
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <memory>
#include <ranges>
#include <typeinfo>
#include <vector>

namespace libjst
{

    template <typename stack_t>
    concept observable_stack = requires
    (std::remove_reference_t<stack_t> & s)
    {
        typename std::remove_cvref_t<stack_t>::reference;

        s.pop();
        { s.top() } -> std::same_as<typename std::remove_cvref_t<stack_t>::reference>;
        s.push(std::declval<typename std::remove_cvref_t<stack_t>::reference>());
    };
    /*!\brief A stack notification registry for observable search algorithms.
     *
     * \details
     *
     * This registry manages a set of libjst::stack_publisher's. During construction a set of observers can register
     * in the registry. Whenever the registry is notified about a stack event (push/pop), it will notify the registered
     * observers. It uses some internal runtime polymorphism which is not exposed observers. Each observer must model
     * the libjst::stack_publisher concept.
     * The attached observers must outlive the registry, otherwise calling a notification event on a dangling reference
     * will lead to undefined behaviour.
     */
    class stack_publisher
    {
    private:

        struct subscriber_base;

        std::vector<std::unique_ptr<subscriber_base>> _subscribers{}; //!< The list of attached observers.

    public:
        /*!\name Constructors, destructor, and assignment
         * \{
         */
        stack_publisher() = default;
        //!\}

        template <observable_stack subscriber_t>
        void subscribe(subscriber_t & subscriber) {
            _subscribers.push_back(std::make_unique<subscriber_impl<subscriber_t>>(subscriber));
        }

        template <observable_stack subscriber_t>
        void unsubscribe(subscriber_t & subscriber) {
            // how do we find this subscriber!
            // subscriber must be able to compare with equality
            std::type_info const & subscriber_type = typeid(subscriber);
            auto it = std::ranges::find_if(_subscribers, [&] (auto && subscriber_ptr) {
                return subscriber_ptr->is_same_type(subscriber_type);
            });
            if (it != std::ranges::end(_subscribers))
                _subscribers.erase(it);
        }

        void notify_push() const
        {
            std::ranges::for_each(_subscribers, [] (auto && subscriber) { subscriber->notify_push(); });
        }

        void notify_pop() const
        {
            std::ranges::for_each(_subscribers, [] (auto && subscriber) { subscriber->notify_pop(); });
        }

    private:
        struct subscriber_base
        {
            subscriber_base() = default; //!< Default.
            virtual ~subscriber_base() = default; //!< Virtual and default.

            //!\brief Notifies the attached observer about a pop event.
            virtual void notify_pop() = 0;
            //!\brief Notifies the attached observer about a push event.
            virtual void notify_push() = 0;
            //!\brief Notifies the attached observer about a push event.
            virtual bool is_same_type(std::type_info const & other_type) const noexcept = 0;
        };

        template <typename subscriber_t>
        struct subscriber_impl final : public subscriber_base
        {
            subscriber_t &_subscriber;

            subscriber_impl(subscriber_t &subscriber) noexcept : _subscriber{subscriber}
            {}

            void notify_pop() noexcept(noexcept(_subscriber.pop())) override
            {
                _subscriber.pop();
            }

            void notify_push() noexcept(noexcept(_subscriber.push(_subscriber.top()))) override
            {
                _subscriber.push(_subscriber.top());
            }

            bool is_same_type(std::type_info const & other_type) const noexcept override
            {
                return typeid(_subscriber) == other_type;
            }
        };
    };

} // namespace libjst
