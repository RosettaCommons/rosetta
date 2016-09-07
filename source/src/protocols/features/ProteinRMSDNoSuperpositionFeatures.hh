// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ProteinRMSDNoSuperpositionFeatures.hh
/// @brief  report ProteinRMSD similarity of structure against supplied reference structure without superposition
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// @author Kyle Barlow (kb@kylebarlow.com)

#ifndef INCLUDED_protocols_features_ProteinRMSDNoSuperpositionFeatures_hh
#define INCLUDED_protocols_features_ProteinRMSDNoSuperpositionFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/ProteinRMSDNoSuperpositionFeatures.fwd.hh>

//External

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace features {

class ProteinRMSDNoSuperpositionFeatures : public FeaturesReporter {
public:

	/// @details If you use this constructor, before applying, make
	/// sure the reference pose is specified, either through the
	/// set_reference_pose method or parse_my_tag.
	ProteinRMSDNoSuperpositionFeatures() {}

	ProteinRMSDNoSuperpositionFeatures(
		core::pose::PoseCOP reference_pose
	);

	ProteinRMSDNoSuperpositionFeatures(
		ProteinRMSDNoSuperpositionFeatures const & ) :
		FeaturesReporter()
	{}

	~ProteinRMSDNoSuperpositionFeatures() override= default;

	/// @brief return string with class name
	std::string
	type_name() const override;

	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const override;

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const override;

private:
	/// @brief generate the protein_rmsd table schema
	void
	write_protein_rmsd_table_schema(
		utility::sql_database::sessionOP db_session) const;

public:
	core::pose::PoseCOP reference_pose() const;
	void reference_pose(core::pose::PoseCOP);


	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & pose) override;

	/// @brief Sets the native pose to the native pose if set, or the starting pose otherwise
	void
	reference_pose_from_options(core::pose::Pose const & pose);

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	) override;

private:
	core::pose::PoseCOP reference_pose_;

};

} // features namespace
} // protocols namespace

#endif // include guard
