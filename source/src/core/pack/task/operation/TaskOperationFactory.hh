// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/TaskOperationFactory.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_TaskOperationFactory_hh
#define INCLUDED_core_pack_task_operation_TaskOperationFactory_hh

// Unit Headers
#include <core/pack/task/operation/TaskOperationFactory.fwd.hh>

// Package Headers
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pack/task/operation/TaskOperationCreator.fwd.hh>
#include <core/pack/task/operation/ResLvlTaskOperationCreator.fwd.hh>
#include <core/pack/task/operation/ResFilterCreator.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.fwd.hh>

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

// singleton class
class TaskOperationFactory
{
public:
	typedef utility::vector1< TaskOperationOP > TaskOperationOPs;
	typedef std::map< std::string, TaskOperationCreatorOP > TaskOperationCreatorMap;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagOP TagOP;
	typedef utility::tag::TagCOP TagCOP;

public:
	static TaskOperationFactory * get_instance();
	void factory_register( TaskOperationCreatorOP );

	/// @brief add a prototype, using its default type name as the map key
	void add_creator( TaskOperationCreatorOP );
	bool has_type( std::string const & ) const;
	/// @brief pass through to child factories
	void add_creator( ResLvlTaskOperationCreatorOP );
	void add_creator( ResFilterCreatorOP );
/// @brief return new TaskOperation by key lookup in task_operation_creator_map_ (new TaskOperation parses Tag if provided)
	TaskOperationOP newTaskOperation(
		std::string const &,
		basic::datacache::DataMap & datamap,
		TagCOP = TagCOP( TagOP( new Tag() ) )
	) const;
/// @brief fills vector with new TaskOperations from nested "TASKOPERATIONS" TagCOP
	void newTaskOperations(	TaskOperationOPs &, basic::datacache::DataMap & datamap, TagCOP ) const;
/// @brief fills vector with new TaskOperations from xml-like tag file
	void newTaskOperations(	TaskOperationOPs &, basic::datacache::DataMap & datamap, std::string const & ) const;

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
	// private constructor/destructor
	TaskOperationFactory();
	virtual ~TaskOperationFactory();
	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static TaskOperationFactory * create_singleton_instance();

private:
#if defined MULTI_THREADED && defined CXX11
	static std::atomic< TaskOperationFactory * > instance_;
#else
	static TaskOperationFactory * instance_;
#endif

	TaskOperationCreatorMap task_operation_creator_map_;
};

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
