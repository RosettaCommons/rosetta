// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ProteinSilentReport.hh
/// @brief  report feature data to database
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_ProteinSilentReport_hh
#define INCLUDED_protocols_features_ProteinSilentReport_hh

// Unit Headers
#include <protocols/features/Report.hh>
#include <protocols/features/ProteinSilentReport.fwd.hh>

//External

// Project Headers
#include <protocols/features/FeaturesReporter.fwd.hh>
#include <protocols/features/ProtocolFeatures.fwd.hh>
#include <protocols/features/BatchFeatures.fwd.hh>
#include <protocols/features/PdbDataFeatures.fwd.hh>
#include <protocols/features/StructureFeatures.fwd.hh>
#include <protocols/features/StructureScoresFeatures.fwd.hh>
#include <protocols/features/ScoreTypeFeatures.fwd.hh>
#include <protocols/features/PoseConformationFeatures.fwd.hh>
#include <protocols/features/PoseCommentsFeatures.fwd.hh>
#include <protocols/features/ProteinResidueConformationFeatures.fwd.hh>
#include <protocols/features/ResidueFeatures.fwd.hh>
#include <protocols/features/ResidueConformationFeatures.fwd.hh>
#include <protocols/features/JobDataFeatures.fwd.hh>
#include <protocols/features/DatabaseFilters.fwd.hh>
#include <protocols/features/FeaturesReporterFactory.fwd.hh>

// Platform Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// C++ Headers
#include <map>
#include <set>

#include <utility/vector1.hh>


namespace protocols {
namespace features {

class ProteinSilentReport : public protocols::features::Report {

public:
	typedef std::map< core::Size, core::pose::PoseOP > StructureMap;

public:
	ProteinSilentReport();

	ProteinSilentReport(ProteinSilentReport const & src);

	~ProteinSilentReport() override;

	core::Size
	version() override;

	bool is_initialized() const;

	void initialize(utility::sql_database::sessionOP db_session);

	void
	apply(
		core::pose::Pose const & pose,
		utility::sql_database::sessionOP db_sesion) override;

	void
	apply(
		core::pose::Pose const & pose,
		utility::sql_database::sessionOP db_sesion,
		std::string const & tag);

	void
	load_pose(
		utility::sql_database::sessionOP db_session,
		StructureID struct_id,
		core::pose::Pose & pose);

	void
	load_pose(
		utility::sql_database::sessionOP db_session,
		StructureID struct_id,
		std::set<core::Size> residue_numbers,
		core::pose::Pose & pose);

	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session ) const;

	void delete_pose_from_tag(
		utility::sql_database::sessionOP db_session,
		std::string const & tag);

	void delete_pose(
		utility::sql_database::sessionOP db_session,
		StructureID const & struct_id);

	core::Size get_protocol_id() const;

	core::Size get_batch_id() const;


private:
	void
	write_protocol_report(
		utility::sql_database::sessionOP db_session
	);

	void write_full_report(
		core::pose::Pose const & pose,
		utility::sql_database::sessionOP db_session,
		std::string const & tag
	);


private:
	void choose_database_filter();

	bool initialized_;
	DatabaseFilterOP database_filter_;

	core::Size protocol_id_;
	core::Size batch_id_;
	StructureMap structure_map_;

	protocols::features::ProtocolFeaturesOP protocol_features_;
	protocols::features::BatchFeaturesOP batch_features_;
	protocols::features::PdbDataFeaturesOP pdb_data_features_;
	protocols::features::StructureFeaturesOP structure_features_;
	protocols::features::StructureScoresFeaturesOP structure_scores_features_;
	protocols::features::ScoreTypeFeaturesOP score_type_features_;
	protocols::features::PoseConformationFeaturesOP pose_conformation_features_;
	protocols::features::PoseCommentsFeaturesOP pose_comments_features_;
	protocols::features::ProteinResidueConformationFeaturesOP protein_residue_conformation_features_;
	protocols::features::ResidueFeaturesOP residue_features_;
	protocols::features::ResidueConformationFeaturesOP residue_conformation_features_;
	protocols::features::JobDataFeaturesOP job_data_features_;
	protocols::features::FeaturesReporterFactory * features_reporter_factory_;
	utility::vector1< protocols::features::FeaturesReporterOP > features_reporters_;
};

} //namespace
} //namespace

#endif //include guard
