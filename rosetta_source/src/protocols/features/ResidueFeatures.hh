// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ResidueFeatures.hh
/// @brief  report residue to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_ResidueFeatures_hh
#define INCLUDED_protocols_features_ResidueFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/ResidueFeatures.fwd.hh>


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

class ResidueFeatures : public protocols::features::FeaturesReporter {
public:
	ResidueFeatures();

	ResidueFeatures(
		core::scoring::ScoreFunctionOP scfxn);

	ResidueFeatures(ResidueFeatures const & src);

	virtual ~ResidueFeatures();

	///@brief return string with class name
	std::string
	type_name() const;

	///@brief return sql statements that setup the right tables
	std::string
	schema() const;

	///@brief return sql statments that add indexes check invariants etc.
	//std::string cleanup() const;

	///@brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		core::Size struct_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_residue_rows(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		core::Size const struct_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_residue_scores_rows(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		core::Size const struct_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_residues_pdb_rows(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		core::Size const struct_id,
		utility::sql_database::sessionOP db_session);

private:

	core::scoring::ScoreFunctionOP scfxn_;

};

} // namespace
} // namespace

#endif // include guard
