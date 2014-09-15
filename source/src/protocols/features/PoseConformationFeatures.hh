// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/PoseConformationFeatures.hh
/// @brief  report comments stored with each pose
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_PoseConformationFeatures_hh
#define INCLUDED_protocols_features_PoseConformationFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/PoseConformationFeatures.fwd.hh>

//External

// Project Headers
#include <core/types.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>
#include <set>

#include <utility/vector1.hh>


namespace protocols{
namespace features{

class PoseConformationFeatures : public protocols::features::FeaturesReporter {
public:
	PoseConformationFeatures(){}

	PoseConformationFeatures(
		PoseConformationFeatures const & ) :
		FeaturesReporter()
	{}


	virtual ~PoseConformationFeatures(){}

	///@brief return string with class name
	std::string
	type_name() const;

	///@brief generate the table schemas and write them to the database
	virtual void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const;

	///@brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	///@brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	void
	delete_record(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	using protocols::features::FeaturesReporter::load_into_pose;

	/// @brief Load information from pose_conformations into pose
	/// @param[in,out] ideal - true if structure is ideal, false if not
	void
	load_into_pose(
		utility::sql_database::sessionOP db_session,
		StructureID struct_id,
		core::pose::Pose & pose,
		bool & ideal
	);

private:

	/// @return true if structure is ideal, false if not
	bool
	load_sequence(
		utility::sql_database::sessionOP db_session,
		StructureID struct_id,
		core::pose::Pose & pose);

	void
	load_fold_tree(
		utility::sql_database::sessionOP db_session,
		StructureID struct_id,
		core::pose::Pose & pose);

	void
	load_jumps(
		utility::sql_database::sessionOP db_session,
		StructureID struct_id,
		core::pose::Pose & pose);

	void
	load_chain_endings(
		utility::sql_database::sessionOP db_session,
		StructureID struct_id,
		core::pose::Pose & pose);


	Size
	report_features_implementation(
		core::pose::Pose const & pose_orig,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	);
};

} // features namespace
} // protocols namespace

#endif // include guard
