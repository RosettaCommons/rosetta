// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/features/RequirementFactory.hh
/// @brief Factory for creating Requirement objects
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifndef INCLUDED_protocols_sewing_requirements_AssemblyRequirementFactory_hh
#define INCLUDED_protocols_sewing_requirements_AssemblyRequirementFactory_hh


// Unit Headers
#include <protocols/sewing/requirements/AssemblyRequirementFactory.fwd.hh>

// Package Headers
#include <protocols/sewing/requirements/AssemblyRequirementCreator.fwd.hh>
#include <protocols/sewing/requirements/AssemblyRequirement.fwd.hh>


// Platform Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/SingletonBase.hh>
// C++ Headers
#include <map>

#include <utility/vector1.hh>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <atomic>
#include <mutex>
#endif
#endif

namespace protocols {
namespace sewing  {
namespace requirements {

/// Create AssemblyRequirements
class AssemblyRequirementFactory: public utility::SingletonBase< AssemblyRequirementFactory > {


public:
	friend class utility::SingletonBase< AssemblyRequirementFactory >;
private:
	// Private constructor to make it singleton managed
	AssemblyRequirementFactory();

	AssemblyRequirementFactory const &
	operator=( AssemblyRequirementFactory const & ); // unimplemented
	//This should be covered by the singleton base
	/*
	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static AssemblyRequirementFactory * create_singleton_instance();
	*/
public:
	// Warning this is not called because of the singleton pattern
	virtual ~AssemblyRequirementFactory();

	//static AssemblyRequirementFactory * get_instance();

	static std::string assembly_requirement_ct_namer( std::string );
	void define_assembly_requirement_subtag( utility::tag::XMLSchemaDefinition & );
	static std::string assembly_requirement_group_name();
	//utility::vector1<std::string> get_all_requirement_names();

	void factory_register(
		AssemblyRequirementCreatorCOP creator
	);

	AssemblyRequirementOP get_requirement(
		std::string const & type_name
	);

	//Again, this should be covered by the singleton base class
	/*
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
	*/
private:
	/*
#if defined MULTI_THREADED
	static std::atomic< AssemblyRequirementFactory * > instance_;
#else
	static AssemblyRequirementFactory * instance_;
#endif
	*/
	typedef std::map< std::string, AssemblyRequirementCreatorCOP > AssemblyRequirementCreatorMap;
	AssemblyRequirementCreatorMap requirement_types_;

};


/// @brief This templated class will register an instance of an
/// RequirementCreator (class T) with the
/// RequirementFactory.  It will ensure that no
/// RequirementCreator is registered twice, and, centralizes this
/// registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class AssemblyRequirementRegistrator : public utility::factory::WidgetRegistrator< AssemblyRequirementFactory, T >
{

public:
	typedef utility::factory::WidgetRegistrator< AssemblyRequirementFactory, T > parent;
	AssemblyRequirementRegistrator() : parent() {}

};


} //namesapce
} //namesapce
} //namesapce

#endif
