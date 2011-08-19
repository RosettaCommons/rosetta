// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/RotamerRecoveryFeatures.hh
/// @brief  report Rotamer Recovery Statistics Scientific Benchmark
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_RotamerRecoveryFeatures_hh
#define INCLUDED_protocols_features_RotamerRecoveryFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/RotamerRecoveryFeatures.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

namespace protocols{
namespace features{

class RotamerRecoveryFeatures : public protocols::features::FeaturesReporter {
public:
	RotamerRecoveryFeatures();

	RotamerRecoveryFeatures(
		core::scoring::ScoreFunctionOP scfxn);

	RotamerRecoveryFeatures(
		RotamerRecoveryFeatures const & src);

	virtual ~RotamerRecoveryFeatures();

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
		utility::vector1< bool > const & relevant_residues,
		core::Size struct_id,
		utility::sql_database::sessionOP db_session);

private:
	core::scoring::ScoreFunctionOP scfxn_;

};

} // namespace
} // namespace

#endif // include guard
