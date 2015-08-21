// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ScoreFunctionFeatures.hh
/// @brief  Add score function parameters to a features database
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_features_ScoreFunctionFeatures_hh
#define INCLUDED_protocols_features_ScoreFunctionFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/ScoreFunctionFeatures.fwd.hh>

//External

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

class ScoreFunctionFeatures : public protocols::features::FeaturesReporter {
public:
	ScoreFunctionFeatures();

	ScoreFunctionFeatures(
		core::scoring::ScoreFunctionOP scfxn,
		std::string const & scfxn_name);

	ScoreFunctionFeatures( ScoreFunctionFeatures const & src );

	virtual ~ScoreFunctionFeatures();

	/// @brief return string with class name
	std::string
	type_name() const;

	/// @brief generate the table schemas and write them to the database
	virtual void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const;

private:
	/// @brief generate the score_function_weights table schema
	void
	write_score_function_weights_table_schema(
		utility::sql_database::sessionOP db_session) const;

public:
	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & /*data*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/);

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	void delete_record(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	);

	void
	insert_score_function_weights_rows(
		core::Size batch_id,
		utility::sql_database::sessionOP db_session
	) const;

private:
	core::scoring::ScoreFunctionOP scfxn_;
	std::string scfxn_name_;

};

} // namespace
} // namespace

#endif // include guard
