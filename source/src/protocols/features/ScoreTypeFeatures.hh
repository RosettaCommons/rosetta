// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ScoreTypeFeatures.hh
/// @brief  Structure scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_ScoreTypeFeatures_hh
#define INCLUDED_protocols_features_ScoreTypeFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/ScoreTypeFeatures.fwd.hh>

//External

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

namespace protocols {
namespace features {

class ScoreTypeFeatures : public protocols::features::FeaturesReporter {
public:
	ScoreTypeFeatures();

	ScoreTypeFeatures(
		core::scoring::ScoreFunctionOP scfxn);

	ScoreTypeFeatures( ScoreTypeFeatures const & src );

	virtual ~ScoreTypeFeatures();

	/// @brief return string with class name
	std::string
	type_name() const;

	/// @brief generate the table schemas and write them to the database
	virtual void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const;

	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	core::Size
	report_features(
		core::pose::Pose const &,
		utility::vector1< bool > const &,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	core::Size
	report_features(
		core::Size const batch_id,
		utility::sql_database::sessionOP db_session);

	void delete_record(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	);

private:
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
