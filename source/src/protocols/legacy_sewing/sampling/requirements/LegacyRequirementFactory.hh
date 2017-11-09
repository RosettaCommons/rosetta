// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/features/LegacyRequirementFactory.hh
/// @brief Factory for creating Requirement objects
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifndef INCLUDED_protocols_legacy_sewing_sampling_requirements_LegacyRequirementFactory_hh
#define INCLUDED_protocols_legacy_sewing_sampling_requirements_LegacyRequirementFactory_hh

// Unit Headers
#include <protocols/legacy_sewing/sampling/requirements/LegacyRequirementFactory.fwd.hh>

// Package Headers
#include <protocols/legacy_sewing/sampling/requirements/LegacyGlobalRequirement.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyIntraSegmentRequirement.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyRequirementCreator.hh>

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
namespace legacy_sewing  {
namespace sampling {
namespace requirements {

/// Create Features Reporters
class LegacyRequirementFactory : public utility::SingletonBase< LegacyRequirementFactory >
{
public:
	friend class utility::SingletonBase< LegacyRequirementFactory >;

private:

	// Private constructor to make it singleton managed
	LegacyRequirementFactory();
	LegacyRequirementFactory(const LegacyRequirementFactory & src) = delete;

	LegacyRequirementFactory const &
	operator=( LegacyRequirementFactory const &  ) = delete;

public:

	// Warning this is not called because of the singleton pattern
	virtual ~LegacyRequirementFactory();

	void factory_register(
		LegacyGlobalRequirementCreatorCOP creator
	);

	void factory_register(
		LegacyIntraSegmentRequirementCreatorCOP creator
	);

	LegacyGlobalRequirementOP get_global_requirement(
		std::string const & type_name
	);

	LegacyIntraSegmentRequirementOP get_intra_segment_requirement(
		std::string const & type_name
	);

	static std::string legacy_global_requirements_ct_namer( std::string const & );
	static std::string legacy_global_requirements_group_name();
	static std::string legacy_intra_segment_requirements_ct_namer( std::string const & );
	static std::string legacy_intra_segment_requirements_group_name();


	void
	define_intra_segment_requirements_subelement( utility::tag::XMLSchemaDefinition & xsd ) const;
	void
	define_global_requirements_subelement( utility::tag::XMLSchemaDefinition & xsd ) const;




private:

	typedef std::map< std::string, LegacyGlobalRequirementCreatorCOP > LegacyGlobalRequirementCreatorMap;
	LegacyGlobalRequirementCreatorMap global_types_;

	typedef std::map< std::string, LegacyIntraSegmentRequirementCreatorCOP > LegacyIntraSegmentRequirementCreatorMap;
	LegacyIntraSegmentRequirementCreatorMap intra_segment_types_;
};


/// @brief This templated class will register an instance of an
/// RequirementCreator (class T) with the
/// LegacyRequirementFactory.  It will ensure that no
/// RequirementCreator is registered twice, and, centralizes this
/// registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class LegacyRequirementRegistrator : public utility::factory::WidgetRegistrator< LegacyRequirementFactory, T >
{

public:
	typedef utility::factory::WidgetRegistrator< LegacyRequirementFactory, T > parent;
	LegacyRequirementRegistrator() : parent() {}

};


} //namesapce
} //namesapce
} //namesapce
} //namespace

#endif
