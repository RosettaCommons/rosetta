// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/StructureFeatures.hh
/// @brief  report Orbital geometry and scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_StructureFeatures_hh
#define INCLUDED_protocols_features_StructureFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/StructureFeatures.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols{
namespace features{

class StructureFeatures : public protocols::features::FeaturesReporter {
public:
	StructureFeatures();

	StructureFeatures( StructureFeatures const & src );

	virtual ~StructureFeatures();

	///@brief return string with class name
	std::string
	type_name() const;

	///@brief return sql statements that setup the right tables
	std::string
	schema() const;

	///@brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	///@brief collect all the feature data for the pose use
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		core::Size protocol_id,
		utility::sql_database::sessionOP db_session
	);

	///@brief collect all the feature data for the pose use
	///This version allows the tag and the input tag to be specificed
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		core::Size protocol_id,
		utility::sql_database::sessionOP db_session,
		std::string const & tag,
		std::string const & input_tag
	);


	void delete_record(
		core::Size struct_id,
		utility::sql_database::sessionOP db_session
	);

	void
	load_into_pose(
		utility::sql_database::sessionOP db_session,
		core::Size struct_id,
		core::pose::Pose & pose);

	void
	load_tag(
		utility::sql_database::sessionOP db_session,
		core::Size struct_id,
		core::pose::Pose & pose);


	core::Size
	get_struct_id(
		utility::sql_database::sessionOP db_session,
		std::string const & tag,
		core::Size const & protocol_id
	);

};

} // namespace
} // namespace

#endif // include guard
