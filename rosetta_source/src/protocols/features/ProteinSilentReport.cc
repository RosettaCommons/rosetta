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
#include <protocols/features/PdbDataFeatures.hh>
#include <protocols/features/PoseCommentsFeatures.hh>
#include <protocols/features/PoseConformationFeatures.hh>
#include <protocols/features/ProteinResidueConformationFeatures.hh>
#include <protocols/features/ResidueConformationFeatures.hh>
#include <protocols/features/JobDataFeatures.hh>
#include <protocols/features/DatabaseFilters.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Platform Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/scoring/ScoreType.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <string>

static basic::Tracer ProteinSilentReportTracer("protocols.features.ProteinSilentReport");

namespace protocols{
namespace features{

using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;
using core::Size;
using core::pose::Pose;
using core::pose::tag_into_pose;
using core::pose::tag_from_pose;

ProteinSilentReport::ProteinSilentReport() :
	initialized_(false),
	database_filter_(NULL),
	protocol_id_(0),
	structure_map_(),
	protocol_features_( new ProtocolFeatures() ),
	pdb_data_features_( new PdbDataFeatures()),
	structure_features_( new StructureFeatures() ),
	structure_scores_features_( new StructureScoresFeatures() ),
	pose_conformation_features_( new PoseConformationFeatures() ),
	pose_comments_features_( new PoseCommentsFeatures() ),
	protein_residue_conformation_features_( new ProteinResidueConformationFeatures() ),
	residue_conformation_features_ (new ResidueConformationFeatures() ),
	job_data_features_ (new JobDataFeatures() )
{
	database_filter_=get_DB_filter_ptr();
}

ProteinSilentReport::ProteinSilentReport(ProteinSilentReport const & src) :
	Report(),
	initialized_(src.initialized_),
	database_filter_(src.database_filter_),
	protocol_id_(src.protocol_id_),
	structure_map_(src.structure_map_),
	protocol_features_(src.protocol_features_),
	pdb_data_features_(src.pdb_data_features_),
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
	pdb_data_features_->write_schema_to_db(db_session);
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
	std::string const & tag
) {

	vector1< bool > relevant_residues(pose.total_residue(), true);

	if (!initialized_){
		write_schema_to_db(db_session);
		if(basic::options::option[basic::options::OptionKeys::out::database_protocol_id].user())
		{
			core::Size protocol_id = basic::options::option[basic::options::OptionKeys::out::database_protocol_id];
			protocol_id_ = protocol_features_->report_features(
				pose, relevant_residues, protocol_id, db_session);
		}else
		{
			protocol_id_ = protocol_features_->report_features(
				pose, relevant_residues, 0, db_session);
		}
		initialized_ = true;

	}

	if(!database_filter_){
		write_full_report(pose,db_session,tag);
		return;
	}

	std::pair<bool, utility::vector1<core::Size> > temp= (*database_filter_)(pose, db_session, protocol_id_);
	bool write_this_pose = temp.first;
	utility::vector1<core::Size> struct_ids_to_delete = temp.second;
	foreach(core::Size struct_id, struct_ids_to_delete){
		delete_pose(db_session,struct_id);
	}
	if(write_this_pose){
		write_full_report(pose,db_session,tag);
	}
}

void
ProteinSilentReport::load_pose(
	sessionOP db_session,
	std::string tag,
	Pose & pose){

	tag_into_pose(pose,tag);

	Size struct_id = structure_features_->get_struct_id(db_session, tag);

	pose_conformation_features_->load_into_pose(db_session, struct_id, pose);
	pdb_data_features_->load_into_pose(db_session,struct_id,pose);
	job_data_features_->load_into_pose(db_session, struct_id, pose);
	pose_comments_features_->load_into_pose(db_session, struct_id, pose);
	protein_residue_conformation_features_->load_into_pose(db_session, struct_id, pose);
	residue_conformation_features_->load_into_pose(db_session,struct_id,pose);

}

bool ProteinSilentReport::is_initialized() const
{
	return initialized_;
}

void
ProteinSilentReport::write_full_report(
	core::pose::Pose const & pose,
	utility::sql_database::sessionOP db_session,
	std::string const & tag) {

	vector1< bool > relevant_residues(pose.total_residue(), true);

	cppdb::transaction transact_guard(*db_session);
	std::string input_tag(protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag());

	Size struct_id = structure_features_->report_features(
		pose, relevant_residues, protocol_id_, db_session, tag,input_tag);

	pose_conformation_features_->report_features(
		pose, relevant_residues, struct_id, db_session);
	pdb_data_features_->report_features(
		pose,relevant_residues,struct_id,db_session);
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

void ProteinSilentReport::delete_pose(utility::sql_database::sessionOP db_session, std::string const & tag)
{
	core::Size struct_id = structure_features_->get_struct_id(db_session,tag);
	delete_pose(db_session,struct_id);
}

void ProteinSilentReport::delete_pose(utility::sql_database::sessionOP db_session, core::Size const & struct_id)
{
	structure_features_->delete_record(struct_id,db_session);
	pose_conformation_features_->delete_record(struct_id,db_session);
	pdb_data_features_->delete_record(struct_id,db_session);
	structure_scores_features_->delete_record(struct_id,db_session);
	pose_comments_features_->delete_record(struct_id,db_session);
	protein_residue_conformation_features_->delete_record(struct_id,db_session);
	residue_conformation_features_->delete_record(struct_id,db_session);
	job_data_features_->delete_record(struct_id,db_session);
}


} //namespace
} //namespace
