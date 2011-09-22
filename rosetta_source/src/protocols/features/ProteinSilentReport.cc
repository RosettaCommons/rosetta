// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ProteinSilentReport.hh
/// @brief  report feature data to database
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ProteinSilentReport.hh>


// Project Headers
#include <protocols/features/FeaturesReporter.hh>

#include <protocols/features/ProtocolFeatures.hh>
#include <protocols/features/StructureFeatures.hh>
#include <protocols/features/StructureScoresFeatures.hh>
#include <protocols/features/PoseCommentsFeatures.hh>
#include <protocols/features/PoseConformationFeatures.hh>
#include <protocols/features/ProteinResidueConformationFeatures.hh>
#include <protocols/features/ResidueConformationFeatures.hh>
#include <protocols/features/JobDataFeatures.hh>

// Platform Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <string>

namespace protocols{
namespace features{


using std::string;
using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;
using core::Size;
using core::pose::Pose;
using core::pose::tag_into_pose;
using core::pose::tag_from_pose;

ProteinSilentReport::ProteinSilentReport() :
	initialized_(false),
	protocol_id_(0),
	structure_map_(),
	protocol_features_( new ProtocolFeatures() ),
	structure_features_( new StructureFeatures() ),
	structure_scores_features_( new StructureScoresFeatures() ),
	pose_conformation_features_( new PoseConformationFeatures() ),
	pose_comments_features_( new PoseCommentsFeatures() ),
	protein_residue_conformation_features_( new ProteinResidueConformationFeatures() ),
	residue_conformation_features_ (new ResidueConformationFeatures() ),
	job_data_features_ (new JobDataFeatures() )
{}

ProteinSilentReport::ProteinSilentReport(ProteinSilentReport const & src) :
	Report(),
	initialized_(src.initialized_),
	protocol_id_(src.protocol_id_),
	structure_map_(src.structure_map_),
	protocol_features_(src.protocol_features_),
	structure_features_(src.structure_features_),
	structure_scores_features_(src.structure_scores_features_),
	pose_conformation_features_(src.pose_conformation_features_),
	pose_comments_features_(src.pose_comments_features_),
	protein_residue_conformation_features_(src.protein_residue_conformation_features_),
	residue_conformation_features_ (src.residue_conformation_features_ ),
	job_data_features_ ( src.job_data_features_)
{}

ProteinSilentReport::~ProteinSilentReport() {}

Size
ProteinSilentReport::version() { return 1; }

void
ProteinSilentReport::write_schema_to_db(
	sessionOP db_session
) const {
	protocol_features_->write_schema_to_db(db_session);
	structure_features_->write_schema_to_db(db_session);
	structure_scores_features_->write_schema_to_db(db_session);
	pose_conformation_features_->write_schema_to_db(db_session);
	pose_comments_features_->write_schema_to_db(db_session);
	protein_residue_conformation_features_->write_schema_to_db(db_session);
	residue_conformation_features_->write_schema_to_db(db_session);
	job_data_features_->write_schema_to_db(db_session);
}


void
ProteinSilentReport::apply(
	Pose const & pose,
	sessionOP db_session){
	apply(pose, db_session, tag_from_pose(pose));
}


void
ProteinSilentReport::apply(
	Pose const & pose,
	sessionOP db_session,
	string const & tag) {

	vector1< bool > relevant_residues(pose.total_residue(), true);

	//db_session->begin();
	cppdb::transaction transact_guard(*db_session);
	if (!initialized_){
		write_schema_to_db(db_session);
		protocol_id_ = protocol_features_->report_features(
		  pose, relevant_residues, 0, db_session);
		initialized_ = true;
	}


	Size struct_id = structure_features_->report_features(
		pose, relevant_residues, protocol_id_, db_session, tag);

	pose_conformation_features_->report_features(
		pose, relevant_residues, struct_id, db_session);
	structure_scores_features_->report_features(
		pose, relevant_residues, struct_id, db_session);
	pose_comments_features_->report_features(
		pose, relevant_residues, struct_id, db_session);
	protein_residue_conformation_features_->report_features(
		pose, relevant_residues, struct_id, db_session);
	residue_conformation_features_->report_features(
		pose, relevant_residues, struct_id,db_session);
	job_data_features_->report_features(
		pose, relevant_residues,struct_id,db_session);

	transact_guard.commit();

}

void
ProteinSilentReport::load_pose(
	sessionOP db_session,
	string tag,
	Pose & pose){

	tag_into_pose(pose,tag);

	Size struct_id = structure_features_->get_struct_id(db_session, tag);

	pose_conformation_features_->load_into_pose(db_session, struct_id, pose);
	pose_comments_features_->load_into_pose(db_session, struct_id, pose);
	protein_residue_conformation_features_->load_into_pose(db_session, struct_id, pose);
	residue_conformation_features_->load_into_pose(db_session,struct_id,pose);

}


} //namespace
} //namespace
