// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ScoreTypeFeatures.hh
/// @brief  Structure scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_ScoreTypeFeatures_hh
#define INCLUDED_protocols_features_ScoreTypeFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/ScoreTypeFeatures.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

namespace protocols{
namespace features{

class ScoreTypeFeatures : public protocols::features::FeaturesReporter {
public:
	ScoreTypeFeatures();

	ScoreTypeFeatures(
		core::scoring::ScoreFunctionOP scfxn);

	ScoreTypeFeatures( ScoreTypeFeatures const & src );

	virtual ~ScoreTypeFeatures();

	///@brief return string with class name
	std::string
	type_name() const;

	///@brief return sql statements that setup the right tables
	std::string
	schema() const;

	///@brief collect all the feature data for the pose
	core::Size
	report_features(
		core::Size protocol_id,
		utility::sql_database::sessionOP db_session
	);

	void delete_record(
		Size struct_id,
		utility::sql_database::sessionOP db_session
	);

	void
	insert_score_type_rows(
		core::Size protocol_id,
		utility::sql_database::sessionOP db_session
	);

private:
	core::scoring::ScoreFunctionOP scfxn_;

};

} // namespace
} // namespace

#endif // include guard
