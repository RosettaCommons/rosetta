// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ProteinResidueConformationFeatures.hh
/// @brief  report idealized torsional DOFs Statistics Scientific Benchmark
/// @author Matthew O'Meara
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_features_ResidueConformationFeatures_hh
#define INCLUDED_protocols_features_ResidueConformationFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/ResidueConformationFeatures.fwd.hh>


// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

namespace protocols{
namespace features{

class ResidueConformationFeatures : public protocols::features::FeaturesReporter {
public:
	ResidueConformationFeatures(){}

	ResidueConformationFeatures(
		ResidueConformationFeatures const & ) :
		FeaturesReporter()
	{}

	virtual ~ResidueConformationFeatures(){}

	///@brief return sql statements that setup the right tables
	std::string
	schema() const;

	///@brief return string with class name
	std::string
	type_name() const;

	///@brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		core::Size struct_id,
		utility::sql_database::sessionOP db_session);

	void
	delete_record(
		core::Size sturct_id,
		utility::sql_database::sessionOP);

	void
	load_into_pose(
		utility::sql_database::sessionOP db_session,
		core::Size struct_id,
		core::pose::Pose & pose);

	void
	load_conformation(
		utility::sql_database::sessionOP db_session,
		core::Size struct_id,
		core::pose::Pose & pose);

private:
	void
	set_coords_for_residue(
		utility::sql_database::sessionOP db_session,
		core::Size struct_id,
		core::Size seqpos,
		core::pose::Pose & pose);


};

} // features
} // protocols

#endif // include guard
