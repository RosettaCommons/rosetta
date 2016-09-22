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


#ifndef INCLUDED_protocols_sewing_sampling_requirements_RequirementFactory_hh
#define INCLUDED_protocols_sewing_sampling_requirements_RequirementFactory_hh

// Unit Headers
#include <protocols/sewing/sampling/requirements/RequirementFactory.fwd.hh>

// Package Headers
#include <protocols/sewing/sampling/requirements/GlobalRequirement.hh>
#include <protocols/sewing/sampling/requirements/IntraSegmentRequirement.hh>
#include <protocols/sewing/sampling/requirements/RequirementCreator.hh>

// Platform Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/factory/WidgetRegistrator.hh>

#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <map>


namespace protocols {
namespace sewing  {
namespace sampling {
namespace requirements {

/// Create Features Reporters
class RequirementFactory : public utility::SingletonBase< RequirementFactory >
{
public:
	friend class utility::SingletonBase< RequirementFactory >;

private:

	// Private constructor to make it singleton managed
	RequirementFactory();
	RequirementFactory(const RequirementFactory & src) = delete;

	RequirementFactory const &
	operator=( RequirementFactory const &  ) = delete;

public:

	// Warning this is not called because of the singleton pattern
	virtual ~RequirementFactory();

	void factory_register(
		GlobalRequirementCreatorCOP creator
	);

	void factory_register(
		IntraSegmentRequirementCreatorCOP creator
	);

	GlobalRequirementOP get_global_requirement(
		std::string const & type_name
	);

	IntraSegmentRequirementOP get_intra_segment_requirement(
		std::string const & type_name
	);

private:

	typedef std::map< std::string, GlobalRequirementCreatorCOP > GlobalRequirementCreatorMap;
	GlobalRequirementCreatorMap global_types_;

	typedef std::map< std::string, IntraSegmentRequirementCreatorCOP > IntraSegmentRequirementCreatorMap;
	IntraSegmentRequirementCreatorMap intra_segment_types_;
};


/// @brief This templated class will register an instance of an
/// RequirementCreator (class T) with the
/// RequirementFactory.  It will ensure that no
/// RequirementCreator is registered twice, and, centralizes this
/// registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class RequirementRegistrator : public utility::factory::WidgetRegistrator< RequirementFactory, T >
{

public:
	typedef utility::factory::WidgetRegistrator< RequirementFactory, T > parent;
	RequirementRegistrator() : parent() {}

};


} //namesapce
} //namesapce
} //namesapce
} //namespace

#endif
