// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rosetta_scripts/PosePropertyReporterFactory.hh
/// @brief Factory for PosePropertyReporters
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifndef INCLUDED_protocols_moves_PosePropertyReporterFactory_hh
#define INCLUDED_protocols_moves_PosePropertyReporterFactory_hh

// Unit Headers
#include <protocols/rosetta_scripts/PosePropertyReporterFactory.fwd.hh>
#include <protocols/rosetta_scripts/PosePropertyReporter.fwd.hh>
#include <protocols/rosetta_scripts/PosePropertyReporter.hh>
#include <protocols/rosetta_scripts/PosePropertyReporterCreator.hh>

// Project Headers
#include <basic/datacache/DataMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/SingletonBase.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// c++ headers
#include <map>
#include <set>

namespace protocols {
namespace rosetta_scripts {

/// @brief This templated class will register an instance of an
/// PosePropertyReporterCreator (class T) with the PosePropertyReporterFactory.  It will ensure
/// that no PosePropertyReporterCreator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class PosePropertyReporterRegistrator : public utility::factory::WidgetRegistrator< PosePropertyReporterFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< PosePropertyReporterFactory, T > parent;
public:
	PosePropertyReporterRegistrator() : parent() {}
};


class PosePropertyReporterFactory: public utility::SingletonBase< PosePropertyReporterFactory >
{
public:
	friend class utility::SingletonBase< PosePropertyReporterFactory >;

	typedef std::map< std::string, PosePropertyReporterCreatorOP > PosePropertyReporterMap;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;
	typedef core::pose::Pose Pose;

public:
	virtual ~PosePropertyReporterFactory();

	void factory_register( PosePropertyReporterCreatorOP creator );

	/// @brief Create a PosePropertyReporter given its identifying string
	PosePropertyReporterOP newPosePropertyReporter( std::string const & );

	/// @brief return new PosePropertyReporter by Tag parsing; the identifying string for the PosePropertyReporter is in the Tag
	PosePropertyReporterOP
	newPosePropertyReporter(
		TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);

	void define_pose_reporter_group( utility::tag::XMLSchemaDefinition & xsd ) const;

	static std::string pose_reporter_group_name();

	static std::string
	complex_type_name_for_pose_reporter( std::string const & reporter_name );

private:
	PosePropertyReporterFactory();

	// Unimplemented -- uncopyable
	PosePropertyReporterFactory( PosePropertyReporterFactory const & ) = delete;
	PosePropertyReporterFactory const & operator = ( PosePropertyReporterFactory const & ) = delete;

private:

	PosePropertyReporterMap reporter_creator_map_;
};

} //namespace rosetta_scripts
} //namespace protocols

#endif
