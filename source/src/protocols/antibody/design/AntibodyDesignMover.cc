// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/design/AntibodyDesignMover.cc
/// @brief Handles the Antibody Design Protocol.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/AntibodyDesignMover.hh>
#include <protocols/antibody/design/AntibodyDesignModeler.hh>
#include <protocols/antibody/database/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/AntibodyGraftDesignMover.hh>
#include <protocols/antibody/design/AntibodySeqDesignMover.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/clusters/util.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/snugdock/SnugDock.hh>
#include <protocols/antibody/snugdock/SnugDockProtocol.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <protocols/antibody/constraints/util.hh>
#include <protocols/antibody/constraints/ParatopeSiteConstraintMover.hh>
#include <protocols/antibody/constraints/ParatopeEpitopeSiteConstraintMover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/jd2/PDBJobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/relax/FastRelax.hh>

//Options
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>


#include <utility/exit.hh>

static thread_local basic::Tracer TR( "protocols.antibody.design.AntibodyDesignMover" );
namespace protocols{
namespace antibody {
namespace design {

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::scoring;
using namespace protocols::antibody;
using namespace core::import_pose;
using namespace protocols::antibody::clusters;

using core::pose::Pose;
using std::string;

AntibodyDesignMover::AntibodyDesignMover() : protocols::moves::Mover(),
		graft_designer_(/* NULL */),
		seq_designer_(/* NULL */),
		modeler_(/* NULL */),
		scorefxn_(/* NULL */),
		design_scorefxn_(/* NULL */),
		ab_info_(/* NULL */)
{
	protocols::moves::Mover::type( "AntibodyDesign" );
	read_cmd_line_options();

}


AntibodyDesignMover::~AntibodyDesignMover(){}

std::string
AntibodyDesignMover::get_name() const {
	return "AntibodyDesignMover";
}

protocols::moves::MoverOP
AntibodyDesignMover::clone() const {
	return protocols::moves::MoverOP( new AntibodyDesignMover(*this) );
}

protocols::moves::MoverOP
AntibodyDesignMover::fresh_instance() const {
	return protocols::moves::MoverOP( new AntibodyDesignMover() );
}
void
AntibodyDesignMover::read_cmd_line_options(){
	set_use_graft_designer(option [OptionKeys::antibody::design::do_graft_design]());
	set_use_sequence_designer(option [OptionKeys::antibody::design::do_sequence_design]());
	set_do_post_graft_design_modeling(option [OptionKeys::antibody::design::do_post_graft_design_modeling]());
	set_do_post_design_modeling(option [OptionKeys::antibody::design::do_post_design_modeling]());
	post_graft_ensemble_output_ = option [OptionKeys::antibody::design::dump_post_graft_designs]();
	interface_distance_ = basic::options::option [basic::options::OptionKeys::antibody::design::interface_dis]();

	//neighbor_detection_dis(basic::options::option [basic::options::OptionKeys::antibody::design::neighbor_dis]());

	if (option [OptionKeys::antibody::design::design_scorefxn].user()){
		set_design_scorefxn(core::scoring::ScoreFunctionFactory::create_score_function(option [OptionKeys::antibody::design::design_scorefxn]()));
	}
	else{
		set_design_scorefxn(get_score_function());
		design_scorefxn_->apply_patch_from_file("antibody_design");
	}

	if (option [OptionKeys::antibody::design::epitope].user()){
		vector1<std::string> epitope_res_strings = option [ OptionKeys::antibody::design::epitope]();
		epitope_residues_ = protocols::antibody::design::get_pdb_numbering_from_string(epitope_res_strings);
	}

	if (basic::options::option [basic::options::OptionKeys::antibody::design::paratope].user()){
		paratope_cdrs_.clear();
		paratope_cdrs_.resize(6, false);
		vector1<std::string> cdrs = basic::options::option [basic::options::OptionKeys::antibody::design::paratope]();
		AntibodyEnumManager manager = AntibodyEnumManager();
		for (core::Size i = 1; i <= cdrs.size(); ++i){
			CDRNameEnum cdr = manager.cdr_name_string_to_enum(cdrs[i]);
			paratope_cdrs_[cdr] = true;
		}
	}
	else{
		paratope_cdrs_.clear();
		paratope_cdrs_.resize(6, true);
	}

}

void
AntibodyDesignMover::setup_scorefxns(){

	if (! scorefxn_){
		scorefxn_ = get_score_function();
		scorefxn_->apply_patch_from_file("antibody_design");
	}
	if (! design_scorefxn_){
		design_scorefxn_ = get_score_function();
		design_scorefxn_->apply_patch_from_file("antibody_design");
	}
}

void
AntibodyDesignMover::setup_constraints(core::pose::Pose& pose){
	for (core::Size i = 1; i<= core::Size(CDRNameEnum_total); ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);

		if (constraints::cdr_has_res_constraints(ab_info_, pose, cdr, "DihedralConstraint")){
			continue;
		} else if (constraints::cdr_has_res_constraints(ab_info_, pose, cdr, "CoordinateConstraint")){
			continue;
		} else {
			constraints::add_harmonic_cluster_cst_or_dihedral_cst(ab_info_, pose, cdr);
		}

	}
}

void
AntibodyDesignMover::set_use_graft_designer(bool setting){
	run_graft_designer_ = setting;
}

void
AntibodyDesignMover::set_use_sequence_designer(bool setting){
	run_sequence_designer_ = setting;
}

void
AntibodyDesignMover::set_do_post_graft_design_modeling(bool setting){
	run_post_graft_modeling_ = setting;
}

void
AntibodyDesignMover::set_do_post_design_modeling(bool setting){
	run_post_design_modeling_ = setting;
}


void
AntibodyDesignMover::set_ab_info(AntibodyInfoOP ab_info){
	ab_info_ = ab_info;
}

void
AntibodyDesignMover::set_scorefxn(ScoreFunctionOP scorefxn){
	scorefxn_ = scorefxn;
}

void
AntibodyDesignMover::set_design_scorefxn(ScoreFunctionOP design_scorefxn){
	design_scorefxn_ = design_scorefxn;
}

void
AntibodyDesignMover::model_post_graft(core::pose::Pose & pose){

	//Look into actual antibody modeling protocols here.
	protocols::moves::MonteCarlo mc = protocols::moves::MonteCarlo(pose, *scorefxn_, .8);

	//Dock + Repack LH interface - Allow Heavy to move relative to L due to H3.
	scorefxn_->show(pose);
	modeler_->ab_dock_chains("L_HA");// Only way this will work currently - LHA pdbs in pose
	modeler_->dock_high_res(pose, true);
	///Add extra high-res step of minimizing the CDRs here instead of doing it after.
	mc.boltzmann(pose);

	//Dock + Repack LH_A interface - Allow Antibody to move relative to antigen.
	modeler_->ab_dock_chains("LH_A");
	modeler_->dock_high_res(pose, true);
	///Add extra high-res step of minimizing the CDRs here instead of doing it after (Think about using CDRsMinPackMin)

	mc.boltzmann(pose);

	//Minimize the cdrs + neighbors to make a tighter fit.  Designed_Relax should allow mutations to fit here.
	//Relaxing the whole structure here may create too tight a fit and negatively impact the SeqDesign phase.
	modeler_->minimize_cdrs(pose);
	mc.boltzmann(pose);
	mc.recover_low(pose);
	TR << "Lowest found: " << (*scorefxn_)(pose) << std::endl;

}

void
AntibodyDesignMover::model_post_design(core::pose::Pose& pose){

	//If no constraints are set - add constraints to pose CDRs.
	core::Real start = scorefxn_->score(pose);

	// Fails at RefineOneCDR.  Debugger is not helping at all.
	SnugDockOP snug( new SnugDock() );
	AntibodyInfoOP updated_ab_info( new AntibodyInfo(pose, AHO_Scheme, North) ); //Should be a reset method in AbInfo..

	snug->set_scorefxn(scorefxn_);
	snug->set_antibody_info(updated_ab_info);//Updated info for movemaps, foldtrees, etc.
	//snug->number_of_high_resolution_cycles(25); //Half as many cycles to save some time as we do a full dualspace relax next.
	//snug->debug();
	snug->apply(pose);
	core::Real post_snugdock = scorefxn_->score(pose);

	//I like relax, and the dualspace protocol is amazing.
	//modeler_->relax_cdrs(pose, true);
	//

	protocols::relax::FastRelaxOP rel( new protocols::relax::FastRelax() );

	ScoreFunctionOP dualspace_scorefxn = scorefxn_->clone(); //May need to be strictly Talaris2013_cart;
	dualspace_scorefxn->set_weight_if_zero(cart_bonded, .5);
	dualspace_scorefxn->set_weight(pro_close, 0);

	rel->set_scorefxn(dualspace_scorefxn);
	rel->dualspace(true);
	rel->max_iter(200);
	rel->minimize_bond_angles(true);
	rel->apply(pose);

	core::Real post_modeling = scorefxn_->score(pose);

	TR <<"start:            " << start <<std::endl;
	TR <<"postSD:        " << post_snugdock << std::endl;
	TR <<"postDsREL: " << post_modeling << std::endl;

}

void
AntibodyDesignMover::output_ensemble(vector1<core::pose::PoseOP> ensemble, core::Size range_start, core::Size range_end, std::string prefix){

	protocols::jd2::JobOP current_job( protocols::jd2::JobDistributor::get_instance()->current_job());
	for (core::Size i = range_start; i <= range_end; ++i){
		//Filter here
		std::string tag = prefix+"_"+utility::to_string(i)+"_";
		TR << "Outputting ensemble " << i << ", "<< tag << std::endl;
		protocols::jd2::JobDistributor::get_instance()->job_outputter()->other_pose(current_job, *(ensemble[i]), tag);
	}
}

void
AntibodyDesignMover::add_cluster_comments_to_pose(core::pose::Pose& pose){

	for (core::SSize i = 1; i <= 6; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);


		if (pose.data().has(core::pose::datacache::CacheableDataType::CDR_CLUSTER_INFO)){
			BasicCDRClusterSet const & cluster_cache = static_cast< BasicCDRClusterSet const & >(pose.data().get(core::pose::datacache::CacheableDataType::CDR_CLUSTER_INFO));
			CDRClusterCOP result = cluster_cache.get_cluster(cdr);
			std::string output = "POSTGRAFT_OR_ORIGINAL_CLUSTER "+ ab_info_->get_cluster_name(result->cluster()) +" "+utility::to_string(result->normalized_distance_in_degrees());
			core::pose::add_comment(pose, "REMARK "+ab_info_->get_CDR_name(cdr), output);
		}

		CDRClusterCOP result = ab_info_->get_CDR_cluster(cdr);
		std::string output = "CLUSTER "+ ab_info_->get_cluster_name(result->cluster()) +" "+utility::to_string(result->normalized_distance_in_degrees());
		core::pose::add_comment(pose, "REMARK "+ab_info_->get_CDR_name(cdr), output);

	}
}

void
AntibodyDesignMover::setup_options_classes() {
	//Overrides Instruction file and defaults.

	cdr_set_options_ = protocols::antibody::design::get_cdr_set_options();
	cdr_graft_design_options_ = protocols::antibody::design::get_graft_design_options();
	cdr_seq_design_options_ = protocols::antibody::design::get_seq_design_options();

	if (option [OptionKeys::antibody::design::design_cdrs].user()){
		utility::vector1<std::string> cdrs = option [OptionKeys::antibody::design::design_cdrs]();
		for (core::Size i = 1; i <= 6; ++i){
			cdr_set_options_[i]->load(false);
			cdr_graft_design_options_[i]->design(false);
			cdr_seq_design_options_[i]->design(false);
		}
		for (core::Size i = 1; i <= cdrs.size(); ++i){
			antibody::CDRNameEnum cdr_enum = ab_info_->get_CDR_name_enum(cdrs[i]);
			cdr_set_options_[cdr_enum]->load(true);
			cdr_graft_design_options_[cdr_enum]->design(true);
			cdr_seq_design_options_[cdr_enum]->design(true);
		}
	}

}

void
AntibodyDesignMover::setup_design_classes() {

	modeler_ = AntibodyDesignModelerOP( new AntibodyDesignModeler(ab_info_) );
	modeler_->set_scorefunction(scorefxn_);
	modeler_->set_cdr_range(CDRNameEnum_start, CDRNameEnum_total, true); //For CDR-CDR interactions

	graft_designer_ = AntibodyGraftDesignMoverOP( new AntibodyGraftDesignMover(ab_info_) );
	graft_designer_->set_scorefunction(scorefxn_);
	graft_designer_->set_epitope_residues(epitope_residues_);
	graft_designer_->set_paratope_cdrs(paratope_cdrs_);

	seq_designer_ = AntibodySeqDesignMoverOP( new AntibodySeqDesignMover(ab_info_) );
	seq_designer_->set_scorefxn(design_scorefxn_);
	seq_designer_->set_epitope_residues(epitope_residues_);
	seq_designer_->set_paratope_cdrs(paratope_cdrs_);

	graft_designer_->set_cdr_set_options(cdr_set_options_);
	graft_designer_->set_graft_design_options(cdr_graft_design_options_);
	graft_designer_->set_seq_design_options(cdr_seq_design_options_);

	seq_designer_->set_seq_designer_options(cdr_seq_design_options_);

}

void
AntibodyDesignMover::setup_epitope_residues(const core::pose::Pose& pose){

	//This is so we setup epitope at the very beginning of the protocol.
	// This needs to also be added to graft designer so that same behavior is used by itself.
	if ( epitope_residues_.size() == 0){
		vector1<bool> epitope = select_epitope_residues(ab_info_, pose, interface_distance_);
		for (core::Size i = 1; i <= pose.total_residue(); ++i){
			if (! epitope[i]) continue;

			PDBNumbering numbering;
			numbering.icode = pose.pdb_info()->icode(i);
			numbering.chain = pose.pdb_info()->chain(i);
			numbering.resnum = pose.pdb_info()->number(i);
			epitope_residues_.push_back(numbering);
		}
	}
}

void
AntibodyDesignMover::reorder_poses(utility::vector1<core::pose::PoseOP>& poses){
    	//Can be refactored to use utility::TopScoreSelector
	//From mc algorithm, you can have multiple poses that are equivalent...

    if (poses.size() <= 1) return;
    utility::vector1<core::pose::PoseOP> sorted_poses;
    sorted_poses.push_back(poses[1]);

    for (core::Size i = 1; i <= poses.size(); ++i){
    	core::Real scU = (*scorefxn_)(*poses[i]);
    	vector1<core::pose::PoseOP>::iterator pose_it = sorted_poses.begin();

    	for (core::Size x = 1; x <= sorted_poses.size(); ++x){
    		core::Real scS = (*scorefxn_)(*sorted_poses[x]);
    		if (scU < scS){
    			sorted_poses.insert(pose_it+x-1, poses[i]);
    			break;
    		}
    		//If the unsorted pose is not lower in energy than any of the sorted poses, add it to the end
    		if (x == sorted_poses.size()){
    			sorted_poses.push_back(poses[i]);
    		}
    	}
    }
    assert(sorted_poses.size() == poses.size());
    poses = sorted_poses;
}


void
AntibodyDesignMover::apply(core::pose::Pose& pose){

	///Setup Objects///


	if (! protocols::antibody::clusters::check_if_pose_renumbered_for_clusters(pose)){
		utility_exit_with_message("PDB must be numbered correctly to identify North CDR clusters.  Please see Antibody Design documentation.");
	}

	if (! ab_info_){
		ab_info_ = AntibodyInfoOP( new AntibodyInfo(pose, AHO_Scheme, North) );
		ab_info_->show(TR);

	}
	if (ab_info_->get_current_AntibodyNumberingScheme()!="AHO_Scheme" && ab_info_->get_current_CDRDefinition() != "North"){
		utility_exit_with_message("Antibody Design Protocol requires AHO_scheme and North definitions");
	}

	setup_scorefxns();
	setup_epitope_residues(pose);
	setup_options_classes();
	setup_design_classes();

	(*scorefxn_)(pose);
	scorefxn_->show(TR, pose);

	//Will go in SHOW.
	TR << std::endl;
	TR << "/// Run Graft Designer: "<< std::boolalpha << run_graft_designer_ <<std::endl;
	TR << "//     ->Post modeling:  "<< std::boolalpha << run_post_graft_modeling_ << std::endl;
	TR << "/////////////////////////" << std::endl;
	TR << "/// Run Seq Designer: " <<std::boolalpha << run_sequence_designer_ << std::endl;
	TR<< "//     -> Post modeling: " << std::boolalpha << run_post_design_modeling_ << std::endl;
	TR << std::endl;

	vector1<core::pose::PoseOP> pose_ensemble;
	vector1<core::pose::PoseOP> final_pose_ensemble;
	if (run_graft_designer_){
		TR <<"Running Graft Design Protocol" <<std::endl;
		graft_designer_->apply(pose);
		pose_ensemble = graft_designer_->get_top_designs();

		for (core::Size i = 1; i <= pose_ensemble.size(); ++i){
			//ab_info_->setup_CDR_clusters(*pose_ensemble[i]);
			//bool removed = pose_ensemble[i]->remove_constraints(); //Should we outright remove all of them?
			//current_constraint_result = protocols::antibody::add_harmonic_cluster_constraints(ab_info_, *pose_ensemble[i]);
			if (run_post_graft_modeling_){
				TR << "Modeling post graft design: ensemble "<< i << std::endl;
				setup_constraints(*pose_ensemble[i]);
				model_post_graft(*pose_ensemble[i]);
			}

			//Use a filter on ensembles generated?  Could have ~15 ensembles. What about the main pose?
			//If it passes filter, add to final_pose_ensemble.
			final_pose_ensemble.push_back(pose_ensemble[i]);
		}
	}
	else{
		//ab_info_->setup_CDR_clusters(pose);
		//current_constraint_result = protocols::antibody::add_harmonic_cluster_constraints(ab_info_, pose);
		final_pose_ensemble.push_back( core::pose::PoseOP( new Pose() ));
		*(final_pose_ensemble[1]) = pose;
	}

	//Optionally output mid-protocol ensembles:
	if (post_graft_ensemble_output_ && run_graft_designer_){
		output_ensemble(final_pose_ensemble, 1, final_pose_ensemble.size(), "mid_ensemble");
	}
	if (run_sequence_designer_){
		TR << "Running sequence designer on pose/ensemble" << std::endl;
		for (core::Size i = 1; i <= final_pose_ensemble.size(); ++i){

			//Quick non-elagent fix so I can run some jobs on the cluster correctly today.
			TR << "Designing ensemble "<< i << "  " << (*scorefxn_)(*final_pose_ensemble[i])<< std::endl;
			seq_designer_->apply(*final_pose_ensemble[i]);
			if (run_post_design_modeling_){
				TR <<"Modeling post protocol" << std::endl;
				setup_constraints(*final_pose_ensemble[i]);
				model_post_design(*final_pose_ensemble[i]);
			}
			TR << "Designed ensemble " << i << std::endl;
			scorefxn_->show(*final_pose_ensemble[i]);

		}
	}

	//Any post-design modeling? Any filters?
	//Reorder pose ensemble before output
	//reorder_poses(final_pose_ensemble); Need debugging - no time to run debugger
	for (core::Size i = 1; i <= final_pose_ensemble.size(); ++i){
		TR << "Pose " << i << ": " << (*scorefxn_)(*final_pose_ensemble[i]) << std::endl;
		ab_info_->setup_CDR_clusters(*final_pose_ensemble[i]);
		add_cluster_comments_to_pose(*final_pose_ensemble[i]);
	}

	//Output final ensembles.
	pose = *final_pose_ensemble[1];
	if (final_pose_ensemble.size() > 1){
		output_ensemble(final_pose_ensemble, 1, final_pose_ensemble.size(), "ensemble");
	}

	if (! option [OptionKeys::out::file::pdb_comments]()){
		TR << "Added pose comments for cluster info.  Use -pdb_comments option to have it output to pdb file." << std::endl;
	}
	TR << "Complete." << std::endl;

}


} //Design
} //Antibody
} //Protocols
