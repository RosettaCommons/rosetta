// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/design/AntibodyDesignMover.cc
/// @brief Handles the Antibody Design Protocol.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/AntibodyDesignMover.hh>
#include <protocols/antibody/design/AntibodyDesignModeler.hh>
#include <protocols/antibody/design/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/AntibodyGraftDesignMover.hh>
#include <protocols/antibody/design/AntibodySeqDesignMover.hh>
#include <protocols/antibody/clusters/util.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/snugdock/SnugDock.hh>
#include <protocols/antibody/snugdock/SnugDockProtocol.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/jd2/PDBJobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

//Options
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/Tracer.hh>


#include <utility/exit.hh>

static basic::Tracer TR("protocols.antibody.design.AntibodyDesignMover");
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
		graft_designer_(NULL),
		seq_designer_(NULL),
		modeler_(NULL),
		scorefxn_(NULL)
{
	ab_info_ = NULL;
	cdr_db_parser_ = new AntibodyDatabaseManager();
	protocols::moves::Mover::type( "AntibodyDesign" );
	read_options();
	set_scorefxn(get_score_function());
}


AntibodyDesignMover::~AntibodyDesignMover(){}
    
std::string
AntibodyDesignMover::get_name() const {
	return "AntibodyDesignMover";
}
    
protocols::moves::MoverOP
AntibodyDesignMover::clone() const {
	return new AntibodyDesignMover(*this);
}

void
AntibodyDesignMover::read_options(){
	set_use_graft_designer(option [OptionKeys::antibody::design::do_graft_design]());
	set_use_sequence_designer(option [OptionKeys::antibody::design::do_sequence_design]());
	set_do_post_graft_design_modeling(option [OptionKeys::antibody::design::do_post_graft_design_modeling]());
	set_do_post_design_modeling(option [OptionKeys::antibody::design::do_post_design_modeling]());
	post_graft_ensemble_output_ = option [OptionKeys::antibody::design::dump_post_graft_designs]();
	
	
	if (option [OptionKeys::antibody::design::design_scorefxn].user()){
		set_design_scorefxn(core::scoring::ScoreFunctionFactory::create_score_function(option [OptionKeys::antibody::design::design_scorefxn]()));
	}
	else{
		set_design_scorefxn(get_score_function());
	}
	
}

void
AntibodyDesignMover::setup_scorefxn(ScoreFunctionOP scorefxn){
	if (scorefxn->get_weight(dihedral_constraint) == 0.0){
		scorefxn->set_weight(dihedral_constraint, 1.0);
	}
	
	if (scorefxn->get_weight(chainbreak) == 0.0){
		scorefxn->set_weight(chainbreak, 100);
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
AntibodyDesignMover::set_graft_designer(AntibodyGraftDesignMoverOP graft_designer){
	graft_designer_ = graft_designer;
}


void
AntibodyDesignMover::set_sequence_designer(AntibodySeqDesignMoverOP seq_designer){
	seq_designer_ = seq_designer;
}

void
AntibodyDesignMover::set_modeler(AntibodyDesignModelerOP modeler){
	modeler_ = modeler;
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
	
	protocols::moves::MonteCarlo mc = protocols::moves::MonteCarlo(pose, *scorefxn_, .8);
	
	//Dock + Repack LH interface - Allow Heavy to move relative to L due to H3.
	scorefxn_->show(pose);
	modeler_->dock_L_H_high_res(pose, true);
	
	mc.boltzmann(pose);
	
	
	
	//Dock + Repack LH_A interface - Allow Antibody to move relative to antigen.
	modeler_->dock_A_LH_high_res(pose, true);
	mc.boltzmann(pose);
	
	
	//Minimize the cdrs + neighbors to make a tighter fit.  Designed_Relax should allow mutations to fit here. 
	modeler_->minimize_cdrs_and_neighbor_sc(pose);
	mc.boltzmann(pose);
	mc.recover_low(pose);
	TR << "Lowest found: " << (*scorefxn_)(pose) << std::endl;
	
}

void
AntibodyDesignMover::model_post_design(core::pose::Pose& pose){
	
	//If no constraints are set - add constraints to pose CDRs.
	
	SnugDockOP snug = new SnugDock();
	snug->set_scorefxn(scorefxn_->clone());
	snug->set_antibody_info(new AntibodyInfo(pose, AHO_Scheme, North));//Updated info for movemaps, foldtrees.
	
	snug->apply(pose);
	
	//I like relax.
	modeler_->relax_cdrs_and_neighbor_sc(pose); 
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
		CDRClusterOP result = ab_info_->get_CDR_cluster(cdr);
		std::string output = "CLUSTER "+ ab_info_->get_cluster_name(result->cluster()) +" "+utility::to_string(result->distance());
		core::pose::add_comment(pose, "REMARK "+ab_info_->get_CDR_name(cdr), output);
		
	}
}

void
AntibodyDesignMover::read_command_line_design_options() {
	//Overrides Instruction file and defaults.
	if (option [OptionKeys::antibody::design::design_cdrs].user()){
		utility::vector1<std::string> cdrs = option [OptionKeys::antibody::design::design_cdrs]();
		for (core::Size i = 1; i <= cdrs.size(); ++i){
			antibody::CDRNameEnum cdr_enum = ab_info_->get_CDR_name_enum(cdrs[i]);
			graft_designer_->set_cdr_range(CDRNameEnum_start, CDRNameEnum_total, false);
			seq_designer_->set_cdr_range(CDRNameEnum_start, CDRNameEnum_total, false);
			
			graft_designer_->set_cdr(cdr_enum, true);
			seq_designer_->set_cdr(cdr_enum, true);
		}
	}
}

void
AntibodyDesignMover::apply(core::pose::Pose& pose){
	
	///Setup Objects///
	

	if (! protocols::antibody::clusters::check_if_pose_renumbered_for_clusters(pose)){
		utility_exit_with_message("PDB must be numbered correctly to identify North CDR clusters.  Please see Antibody Design documentation.");
	}
	
	if (! ab_info_){
		ab_info_ = new AntibodyInfo(pose, AHO_Scheme, North);
		ab_info_->show(TR);

	}
	if (ab_info_->get_current_AntibodyNumberingScheme()!="AHO_Scheme" && ab_info_->get_current_CDRDefinition() != "North"){
		utility_exit_with_message("Antibody Design Protocol requires AHO_scheme and North definitions");
	}
	
	setup_scorefxn(scorefxn_);
	setup_scorefxn(design_scorefxn_);
	
	if (!modeler_ ){
		modeler_ = new AntibodyDesignModeler(ab_info_);
		modeler_->set_scorefunction(scorefxn_);
		modeler_->set_cdr_range(CDRNameEnum_start, CDRNameEnum_total, true); //For CDR-CDR interactions
	}
	
	if (! graft_designer_ ){
		graft_designer_ = new AntibodyGraftDesignMover(ab_info_);
		graft_designer_->set_scorefunction(scorefxn_);
	}
	
	if (! seq_designer_ ){
		seq_designer_ = new AntibodySeqDesignMover(ab_info_);
		seq_designer_->set_scorefxn(design_scorefxn_);
	}
    
	read_command_line_design_options();
	
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
		TR <<"Running Graft Design" <<std::endl;
		graft_designer_->apply(pose);
		pose_ensemble = graft_designer_->get_top_designs();
		
		for (core::Size i = 1; i <= pose_ensemble.size(); ++i){
			//ab_info_->setup_CDR_clusters(*pose_ensemble[i]);
			//bool removed = pose_ensemble[i]->remove_constraints(); //Should we outright remove all of them?
			//current_constraint_result = protocols::antibody::add_harmonic_cluster_constraints(ab_info_, *pose_ensemble[i]);
			if (run_post_graft_modeling_){
				TR << "Modeling post graft design: ensemble "<< i << std::endl;
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
		final_pose_ensemble.push_back( new Pose());
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
				model_post_design(*final_pose_ensemble[i]);
			}
			TR << "Designed ensemble " << i << std::endl;
			scorefxn_->show(*final_pose_ensemble[i]);
			
		}
	}
	
	//Any post-design modeling? Any filters?
	//Reorder pose ensemble before output
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
