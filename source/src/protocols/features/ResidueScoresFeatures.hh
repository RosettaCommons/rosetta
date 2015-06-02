// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ResidueScoresFeatures.hh
/// @brief  report residue scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara (mattjomeara@gmail.com

#ifndef INCLUDED_protocols_features_ResidueScoresFeatures_hh
#define INCLUDED_protocols_features_ResidueScoresFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/ResidueScoresFeatures.fwd.hh>

//External

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <string>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>


namespace protocols{
namespace features{

class ResidueScoresFeatures : public protocols::features::FeaturesReporter {
public:
	ResidueScoresFeatures();

	ResidueScoresFeatures(
		core::scoring::ScoreFunctionOP scfxn);

	ResidueScoresFeatures(ResidueScoresFeatures const & src);

	virtual ~ResidueScoresFeatures();

	/// @brief return string with class name
	std::string
	type_name() const;

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

private:
	/// @brief generate the residue_scores_1b table schema
	void
	write_residue_scores_1b_table_schema(
		utility::sql_database::sessionOP db_session) const;

	///@brief 2b schema helper
	void
	write_residue_scores_2b_table_schema_helper(
		std::string name,
		utility::sql_database::sessionOP db_session) const;
	
	///@brief generate the residue_scores_2b table schema
	void
	write_residue_scores_2b_table_schema(
		utility::sql_database::sessionOP db_session) const;

	/// @brief generate the residue_scores_2b table schema
	void
	write_residue_scores_lr_2b_table_schema(
		utility::sql_database::sessionOP db_session) const;

public:
	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
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

private:

	void
	insert_residue_scores_rows(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID const struct_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_one_body_residue_score_rows(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		Size const batch_id,
		StructureID const struct_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_two_body_residue_score_rows(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		Size const batch_id,
		StructureID const struct_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_two_body_long_range_residue_score_rows(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		Size const batch_id,
		StructureID const struct_id,
		utility::sql_database::sessionOP db_session);

private:

	core::scoring::ScoreFunctionOP scfxn_;

};

} // namespace
} // namespace

#endif // include guard
