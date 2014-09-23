// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/ResFilterFactory.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_ResFilterFactory_hh
#define INCLUDED_core_pack_task_operation_ResFilterFactory_hh

// Unit Headers
#include <core/pack/task/operation/ResFilterFactory.fwd.hh>

// Package Headers
#include <core/pack/task/operation/ResFilter.fwd.hh>
#include <core/pack/task/operation/ResFilterCreator.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>

// c++ headers
#include <string>
#include <map>

#include <utility/vector0.hh>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <atomic>
#include <mutex>
#endif
#endif

namespace core {
namespace pack {
namespace task {
namespace operation {

class ResFilterFactory // singleton; no need to derive from RefCount : public utility::pointer::ReferenceCount
{
public:
	//typedef utility::pointer::ReferenceCount parent;
	typedef std::map< std::string, ResFilterCreatorOP > ResFilterCreatorMap;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;

public:
	static ResFilterFactory * get_instance();
	void factory_register( ResFilterCreatorOP );

	///@brief add a prototype, using its default type name as the map key
	void add_creator( ResFilterCreatorOP );
	bool has_type( std::string const & ) const;

///@brief return new ResFilter by key lookup in filter_map_ (new ResFilter parses Tag if provided)
	ResFilterOP newResFilter( std::string const &, TagCOP = new Tag() ) const;

#ifdef MULTI_THREADED
#ifdef CXX11
public:

	/// @brief This public method is meant to be used only by the
	/// utility::thread::safely_create_singleton function and not meant
	/// for any other purpose.  Do not use.
	static std::mutex & singleton_mutex();

private:
	static std::mutex singleton_mutex_;
#endif
#endif

private:
	ResFilterFactory();
	virtual ~ResFilterFactory();
	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static ResFilterFactory * create_singleton_instance();

private:
#if defined MULTI_THREADED && defined CXX11
	static std::atomic< ResFilterFactory * > instance_;
#else
	static ResFilterFactory * instance_;
#endif

	ResFilterCreatorMap filter_creator_map_;

};

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
