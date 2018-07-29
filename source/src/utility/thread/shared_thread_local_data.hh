// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/thread/shared_thread_local_data.hh
/// @brief  SharedThreadLocalData class, - container for data that needs to instantiated as single, shared instance per-thread
/// @author Sergey Lyskov (sergey.lyskov@jhu.edu)

#ifndef INCLUDED_utility_thread_shared_thread_local_data_hh
#define INCLUDED_utility_thread_shared_thread_local_data_hh

#include <utility/pointer/owning_ptr.hh>

#include <thread>

namespace utility {
namespace thread {

/// This template class is intended to provide access to instance of T that is unique to your thread
/// Notes:
/// - underlying object will be create on-demand. You should not try to store it or share it (thats why we expose only reference to object and not SP)
///
/// - class itself is not thread-safe!
///
/// - you can safely move instance of SharedThreadLocalData between threads, such move will result in creation of underlying object next time it is accessed.
///
/// - So if your use-pattern involved of constantly moving-between-threads it will create two types of overheads:
///
///   - if thread to which you moved your object ALREADY HAVE other instances of SharedThreadLocalData for that type and such instance already hold a valid underlying data (on which `get` was called):
///     overhead will be minimal: one lock + map lookup during the first access and no overhead (except the cost of `std::this_thread::get_id()`) on sub-sequential accesses
///
///   - if thread to which you moved your object DOES NOT HAVE other instances of SharedThreadLocalData then overhead will be:
///     cost of underlying object creation during first access and no overhead (except the cost of `std::this_thread::get_id()`) on sub-sequential accesses
///
///
///  WARNING WARNING WARNING:
///  In the unlikely event that you actually do need to use this class please make to contact Sergey Lyskov (sergey.lyskov@jhu.edu) or Andrew Leaver-Fay (aleaverfay@gmail.com)
///
template <typename T>
class SharedThreadLocalData
{
public:

	T & get();

private:
	std::thread::id thread;  // store thread-id in which object was created
	utility::pointer::shared_ptr<T> data;

private:
	static utility::pointer::shared_ptr<T> get_instance(std::thread::id);
};



} // thread
} // utility


#endif // INCLUDED_utility_thread_shared_thread_local_data_hh
