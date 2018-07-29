// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/thread/shared_thread_local_data.hh
/// @brief  Implementation of SharedThreadLocalData class, - container for data that needs to instantiated as single, shared instance per-thread
/// @author Sergey Lyskov (sergey.lyskov@jhu.edu)


#ifndef INCLUDED_utility_thread_shared_thread_local_data_impl_hh
#define INCLUDED_utility_thread_shared_thread_local_data_impl_hh

#include <utility/thread/shared_thread_local_data.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

#include <mutex>
#include <map>

namespace utility {
namespace thread {

template <typename T>
T & SharedThreadLocalData<T>::get()
{
	auto current_thread_id = std::this_thread::get_id();

	if ( thread != current_thread_id ) data = get_instance(current_thread_id);

	return *data;
}


template <typename T>
utility::pointer::shared_ptr<T> SharedThreadLocalData<T>::get_instance(std::thread::id thread)
{
	static std::mutex mutex;

	static std::map<std::thread::id, utility::pointer::weak_ptr<T> > objects;

	std::lock_guard<std::mutex> guard(mutex);

	utility::pointer::weak_ptr<T> & o_wp = objects[thread];

	if ( auto o_sp = o_wp.lock() ) return o_sp;
	else {
		// we DO NOT use `make_shared` here on purpose: we want to split allocated memory segments for shared_ptr and T, so they belong to two distict allocations
		// this will allow memory allocated for T to be feed when only weak_ptr remain
		o_sp = utility::pointer::shared_ptr<T>( new T() );

		o_wp = o_sp;

		return o_sp;
	}
}



} // thread
} // utility


#endif // INCLUDED_utility_thread_shared_thread_local_data_impl_hh
