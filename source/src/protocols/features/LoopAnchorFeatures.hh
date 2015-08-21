// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/LoopAnchorFeatures.hh
/// @brief  report comments stored with each pose
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_LoopAnchorFeatures_hh
#define INCLUDED_protocols_features_LoopAnchorFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/LoopAnchorFeatures.fwd.hh>

// Project Headers
#include <core/conformation/Residue.fwd.hh>

//External
#include <cppdb/frontend.h>

// Utility Headers
#include <numeric/HomogeneousTransform.fwd.hh>

namespace protocols {
namespace features {

class LoopAnchorFeatures : public FeaturesReporter {
public:
	LoopAnchorFeatures();

	LoopAnchorFeatures(LoopAnchorFeatures const & );

	virtual ~LoopAnchorFeatures();

	/// @brief return string with class name
	std::string
	type_name() const;

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

private:
	/// @brief generate the loop_anchors table schema
	void
	write_loop_anchors_table_schema(
		utility::sql_database::sessionOP db_session) const;

	/// @brief generate the loop_anchor_transforms table schema
	void
	write_loop_anchor_transforms_table_schema(
		utility::sql_database::sessionOP db_session) const;

	/// @brief generate the loop_anchor_transforms_three_res table schema
	void
	write_loop_anchor_transforms_three_res_table_schema(
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
		utility::vector1< bool > const & /*relevant_residues*/,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	void
	set_use_relevant_residues_as_loop_length( bool const use_relevant_residues_as_loop_length );

	void
	set_use_single_residue_to_define_anchor_transfrom( bool const use_single_residue_to_define_anchor_transfrom );

private:
	core::Size min_loop_length( utility::vector1< bool > const & relevant_residue ) const;
	core::Size max_loop_length( utility::vector1< bool > const & relevant_residue ) const;
	core::Size determine_correct_length( utility::vector1< bool > const & relevant_residue, Size default_length ) const;

	numeric::HomogeneousTransform<core::Real> const
	frame_for_residue(core::conformation::Residue const & residue) const;

	numeric::HomogeneousTransform<core::Real> const
	compute_anchor_transform(
		core::pose::Pose const & pose,
		core::Size const residue_begin,
		core::Size const residue_end) const;

	void
	compute_transform_and_write_to_db(
		StructureID struct_id,
		core::Size begin,
		core::Size end,
		core::pose::Pose const & pose,
		cppdb::statement & stmt) const;

private:

	bool use_relevant_residues_as_loop_length_;
	bool use_single_residue_to_define_anchor_transfrom_;
	core::Size min_loop_length_;
	core::Size max_loop_length_;

};

} // features namespace
} // protocols namespace

#endif //INCLUDED_protocols_features_LoopAnchorFeatures_hh
