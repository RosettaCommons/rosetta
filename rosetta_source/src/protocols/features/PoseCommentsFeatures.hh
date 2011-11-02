// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/PoseCommentsFeatures.hh
/// @brief  report comments stored with each pose
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_PoseCommentsFeatures_hh
#define INCLUDED_protocols_features_PoseCommentsFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/PoseCommentsFeatures.fwd.hh>

// Project Headers
#include <core/types.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.fwd.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols{
namespace features{

class PoseCommentsFeatures : public FeaturesReporter {
public:
	PoseCommentsFeatures(){}

	PoseCommentsFeatures(
		PoseCommentsFeatures const & ) :
		FeaturesReporter()
	{}

	virtual ~PoseCommentsFeatures(){}

	///@brief return string with class name
	std::string
	type_name() const;

	///@brief return sql statements that setup the right tables
	std::string
	schema() const;

	///@brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & /*relevant_residues*/,
		core::Size struct_id,
		utility::sql_database::sessionOP db_session);

	void delete_record(
		core::Size struct_id,
		utility::sql_database::sessionOP db_session);

	void
	load_into_pose(
		utility::sql_database::sessionOP db_session,
		core::Size struct_id,
		core::pose::Pose & pose);

};

} // features namespace
} // protocols namespace

#endif // include guard
