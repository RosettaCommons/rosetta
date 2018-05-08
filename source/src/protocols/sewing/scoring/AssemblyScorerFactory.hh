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


#ifndef INCLUDED_protocols_sewing_scoring_AssemblyScorerFactory_hh
#define INCLUDED_protocols_sewing_scoring_AssemblyScorerFactory_hh


// Unit Headers
#include <protocols/sewing/scoring/AssemblyScorerFactory.fwd.hh>

// Package Headers
#include <protocols/sewing/scoring/AssemblyScorerCreator.fwd.hh>
#include <protocols/sewing/scoring/AssemblyScorer.fwd.hh>


// Platform Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/factory/WidgetRegistrator.hh>
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
namespace scoring {

/// Create AssemblyScorers
class AssemblyScorerFactory: public utility::SingletonBase< AssemblyScorerFactory > {


public:
	friend class utility::SingletonBase< AssemblyScorerFactory >;
private:
	// Private constructor to make it singleton managed
	AssemblyScorerFactory();

	AssemblyScorerFactory const &
	operator=( AssemblyScorerFactory const & ); // unimplemented
	/*
	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static AssemblyScorerFactory * create_singleton_instance();
	*/
public:
	// Warning this is not called because of the singleton pattern
	virtual ~AssemblyScorerFactory();

	//static AssemblyScorerFactory * get_instance();

	void factory_register(
		AssemblyScorerCreatorCOP creator
	);

	AssemblyScorerOP get_assembly_scorer(
		std::string const & type_name
	);
	static std::string assembly_scorer_ct_namer( std::string tag_name );
	void define_assembly_scorer_subtag( utility::tag::XMLSchemaDefinition & );
	static std::string assembly_scorer_group_name();

	//utility::vector1<std::string> get_all_assembly_scorer_names();
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
#if defined MULTI_THREADED && defined CXX11
	static std::atomic< AssemblyScorerFactory * > instance_;
#else
	static AssemblyScorerFactory * instance_;
#endif
	*/

	typedef std::map< std::string, AssemblyScorerCreatorCOP > AssemblyScorerCreatorMap;
	AssemblyScorerCreatorMap assembly_scorer_types_;

};


/// @brief This templated class will register an instance of an
/// RequirementCreator (class T) with the
/// RequirementFactory.  It will ensure that no
/// RequirementCreator is registered twice, and, centralizes this
/// registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class AssemblyScorerRegistrator : public utility::factory::WidgetRegistrator< AssemblyScorerFactory, T >
{

public:
	typedef utility::factory::WidgetRegistrator< AssemblyScorerFactory, T > parent;
	AssemblyScorerRegistrator() : parent() {}

};


} //namesapce
} //namesapce
} //namesapce

#endif
