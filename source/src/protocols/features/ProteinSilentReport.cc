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

// Unit Headers
#include <protocols/features/ProteinSilentReport.hh>

//External

// Project Headers
#include <protocols/features/FeaturesReporter.hh>

#include <protocols/features/ProtocolFeatures.hh>
#include <protocols/features/BatchFeatures.hh>
#include <protocols/features/StructureFeatures.hh>
#include <protocols/features/ScoreTypeFeatures.hh>
#include <protocols/features/StructureScoresFeatures.hh>
#include <protocols/features/PdbDataFeatures.hh>
#include <protocols/features/PoseCommentsFeatures.hh>
#include <protocols/features/PoseConformationFeatures.hh>
#include <protocols/features/ProteinResidueConformationFeatures.hh>
#include <protocols/features/ResidueFeatures.hh>
#include <protocols/features/ResidueConformationFeatures.hh>
#include <protocols/features/JobDataFeatures.hh>
#include <protocols/features/DatabaseFilters.hh>
#include <protocols/features/FeaturesReporterFactory.hh>
#include <protocols/features/util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <basic/options/option.hh>
#include <basic/database/sql_utils.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Platform Headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/util.hh>

#include <core/scoring/ScoreType.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// Boost Headers
#include <boost/foreach.hpp>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <string>

#include <protocols/jd2/Job.hh>

namespace protocols {
namespace features {

static THREAD_LOCAL basic::Tracer TR( "protocols.features.ProteinSilentReport" );

using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;
using core::Size;
using core::pose::Pose;
using core::pose::tag_into_pose;
using core::pose::tag_from_pose;

ProteinSilentReport::ProteinSilentReport() :
	initialized_(false),
	database_filter_(/* NULL */),
	protocol_id_(0),
	batch_id_(0),
	structure_map_(),
	protocol_features_( protocols::features::ProtocolFeaturesOP( new ProtocolFeatures() ) ),
	batch_features_( protocols::features::BatchFeaturesOP( new BatchFeatures() ) ),
	pdb_data_features_( protocols::features::PdbDataFeaturesOP( new PdbDataFeatures() )),
	structure_features_( protocols::features::StructureFeaturesOP( new StructureFeatures() ) ),
	structure_scores_features_( protocols::features::StructureScoresFeaturesOP( new StructureScoresFeatures() ) ),
	score_type_features_(protocols::features::ScoreTypeFeaturesOP( new ScoreTypeFeatures() )),
	pose_conformation_features_( protocols::features::PoseConformationFeaturesOP( new PoseConformationFeatures() ) ),
	pose_comments_features_( protocols::features::PoseCommentsFeaturesOP( new PoseCommentsFeatures() ) ),
	protein_residue_conformation_features_( protocols::features::ProteinResidueConformationFeaturesOP( new ProteinResidueConformationFeatures() ) ),
	residue_features_ (protocols::features::ResidueFeaturesOP( new ResidueFeatures() ) ),
	residue_conformation_features_ (protocols::features::ResidueConformationFeaturesOP( new ResidueConformationFeatures() ) ),
	job_data_features_ (protocols::features::JobDataFeaturesOP( new JobDataFeatures() ) )
{
	database_filter_=get_DB_filter_ptr();

	if ( basic::options::option[basic::options::OptionKeys::out::database_protocol_id].user() ) {
		protocol_id_ = basic::options::option[basic::options::OptionKeys::out::database_protocol_id];
	}
	features_reporters_.push_back(structure_features_);
	features_reporters_.push_back(pdb_data_features_);
	features_reporters_.push_back(score_type_features_);
	features_reporters_.push_back(pose_conformation_features_);
	features_reporters_.push_back(pose_comments_features_);
	features_reporters_.push_back(protein_residue_conformation_features_);
	features_reporters_.push_back(residue_features_);
	features_reporters_.push_back(residue_conformation_features_);
	features_reporters_.push_back(job_data_features_);
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
	score_type_features_(src.score_type_features_),
	pose_conformation_features_(src.pose_conformation_features_),
	pose_comments_features_(src.pose_comments_features_),
	protein_residue_conformation_features_(src.protein_residue_conformation_features_),
	residue_features_(src.residue_features_),
	residue_conformation_features_ (src.residue_conformation_features_ ),
	job_data_features_ ( src.job_data_features_)
{}

ProteinSilentReport::~ProteinSilentReport() = default;

Size
ProteinSilentReport::version() { return 1; }

void
ProteinSilentReport::write_schema_to_db(
	sessionOP db_session
) const {
	protocol_features_->write_schema_to_db(db_session, protocol_id_);
	batch_features_->write_schema_to_db(db_session, batch_id_);
	structure_features_->write_schema_to_db(db_session);
	score_type_features_->write_schema_to_db(db_session);
	structure_scores_features_->write_schema_to_db(db_session);
	pdb_data_features_->write_schema_to_db(db_session);
	pose_conformation_features_->write_schema_to_db(db_session);
	pose_comments_features_->write_schema_to_db(db_session);
	residue_features_->write_schema_to_db(db_session);
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

	if ( pose.energies().energies_updated() ) {
		features_reporters_.push_back(structure_scores_features_);
	}

	vector1< bool > relevant_residues(pose.total_residue(), true);

	initialize(db_session);

	std::string input_tag(protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag());
	structure_features_->mark_structure_as_sampled(batch_id_,tag,input_tag,db_session);

	if ( !database_filter_ ) {
		write_full_report(pose,db_session,tag);
		return;
	}

	std::pair<bool, utility::vector1<StructureID> > temp= (*database_filter_)(pose, db_session, protocol_id_);
	bool write_this_pose = temp.first;
	utility::vector1<StructureID> struct_ids_to_delete = temp.second;
	BOOST_FOREACH ( StructureID struct_id, struct_ids_to_delete ) {
		delete_pose(db_session,struct_id);
	}
	if ( write_this_pose ) {
		write_full_report(pose,db_session,tag);
	}
}

void
ProteinSilentReport::load_pose(
	sessionOP db_session,
	StructureID struct_id,
	Pose & pose){

	bool ideal = true; // Set by load_into_pose, and then used by protein_residue_conformation_features to determine
	// if backbone torsions should be loaded into pose
	pose_conformation_features_->load_into_pose(db_session, struct_id, pose, ideal);
	pdb_data_features_->load_into_pose(db_session,struct_id,pose);
	job_data_features_->load_into_pose(db_session, struct_id, pose);
	pose_comments_features_->load_into_pose(db_session, struct_id, pose);
	protein_residue_conformation_features_->load_into_pose(db_session, struct_id, pose, ideal);
	residue_conformation_features_->load_into_pose(db_session,struct_id,pose);
}

void
ProteinSilentReport::load_pose(
	sessionOP db_session,
	StructureID struct_id,
	std::set<core::Size> residue_numbers,
	Pose & pose){

	//first load in the entire pose, then delete unspecified residues. This will be slower than
	//it needs to be but should result in the most sensible state of the fold/atom trees.
	load_pose(db_session, struct_id, pose);

	int total_res=pose.total_residue();
	int num_removed_residues=0;
	for ( int i=1; i<=total_res; ++i ) {
		//if the residue wasn't specified, delete it
		if ( residue_numbers.find(i)==residue_numbers.end() ) {
			pose.conformation().delete_residue_slow(i-num_removed_residues);
			++num_removed_residues;
		}
	}
}

bool ProteinSilentReport::is_initialized() const
{
	return initialized_;
}

void ProteinSilentReport::initialize(sessionOP db_session)
{
	if ( !initialized_ ) {
		write_schema_to_db(db_session);
		write_protocol_report(db_session);
		initialized_ = true;
	}
}

void
ProteinSilentReport::write_protocol_report(
	utility::sql_database::sessionOP db_session
){
	//initialize protocol and batch id
	std::pair<Size, Size> ids = get_protocol_and_batch_id("db_job_outputter", "", features_reporters_, db_session);
	protocol_id_=ids.first;
	batch_id_=ids.second;

	db_session->begin_transaction();
	score_type_features_->report_features(batch_id_, db_session);
	db_session->commit_transaction();
}
void
ProteinSilentReport::write_full_report(
	core::pose::Pose const & pose,
	utility::sql_database::sessionOP db_session,
	std::string const & tag
) {

	vector1< bool > relevant_residues(pose.total_residue(), true);

	db_session->begin_transaction();
	std::string input_tag(protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag());

	StructureID struct_id = structure_features_->report_features(
		batch_id_, db_session, tag, input_tag);

	TR.Info << "Beginning report, struct_id: " << struct_id << std::endl;

	pose_conformation_features_->report_features(
		pose, relevant_residues, struct_id, db_session);

	pdb_data_features_->report_features(
		pose,relevant_residues,struct_id,db_session);

	if ( pose.energies().energies_updated() ) {
		structure_scores_features_->report_features(
			pose, relevant_residues, struct_id, db_session);
	}

	pose_comments_features_->report_features(
		pose, relevant_residues, struct_id, db_session);
	residue_features_->report_features(
		pose, relevant_residues, struct_id, db_session);
	protein_residue_conformation_features_->report_features(
		pose, relevant_residues, struct_id, db_session);
	residue_conformation_features_->report_features(
		pose, relevant_residues, struct_id, db_session);
	job_data_features_->report_features(
		pose, relevant_residues, struct_id, db_session);

	db_session->commit_transaction();
}

void ProteinSilentReport::delete_pose_from_tag(utility::sql_database::sessionOP db_session, std::string const & tag)
{
	StructureID struct_id = structure_features_->get_struct_id(db_session,tag,protocol_id_);
	delete_pose(db_session,struct_id);
}

void ProteinSilentReport::delete_pose(utility::sql_database::sessionOP db_session, StructureID const & struct_id)
{

	job_data_features_->delete_record(struct_id,db_session);
	residue_conformation_features_->delete_record(struct_id,db_session);
	protein_residue_conformation_features_->delete_record(struct_id,db_session);
	residue_features_->delete_record(struct_id,db_session);
	pose_comments_features_->delete_record(struct_id,db_session);
	structure_scores_features_->delete_record(struct_id,db_session);
	pdb_data_features_->delete_record(struct_id,db_session);
	pose_conformation_features_->delete_record(struct_id,db_session);
	structure_features_->delete_record(struct_id,db_session);
}

core::Size ProteinSilentReport::get_protocol_id() const
{
	return protocol_id_;
}

core::Size ProteinSilentReport::get_batch_id() const
{
	return batch_id_;
}


} //namespace
} //namespace
