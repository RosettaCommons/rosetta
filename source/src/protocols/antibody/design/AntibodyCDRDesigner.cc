// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody_design/AntibodyCDRDesigner.cc
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/AntibodyCDRDesigner.hh>
#include <protocols/antibody/design/ConservativeDesignOperation.hh>
#include <protocols/antibody/design/ResidueProbDesignOperation.hh>
#include <protocols/antibody/design/DesignInstructionsParser.hh>
#include <protocols/antibody/design/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/clusters/util.hh>
#include <protocols/antibody/util.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/util.hh>

#include <protocols/relax/FastRelax.hh>
#include <protocols/toolbox/task_operations/LimitAromaChi2Operation.hh>
#include <protocols/toolbox/task_operations/RestrictToLoopsAndNeighbors.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/flxbb/FlxbbDesign.hh>
#include <protocols/flxbb/DesignTask.hh>
#include <protocols/moves/MonteCarlo.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
//#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.antibody.design.AntibodyCDRDesigner");

namespace protocols {
namespace antibody {
namespace design {
	using namespace protocols::antibody;
	using namespace protocols::antibody::clusters;
	
	using namespace core::pack::task::operation;
	using namespace protocols::toolbox::task_operations;

AntibodyCDRDesigner::AntibodyCDRDesigner(AntibodyInfoOP ab_info){
	ab_info_ = ab_info;
	if (ab_info_->get_current_AntibodyNumberingScheme()!="AHO_Scheme" && ab_info_->get_current_CDRDefinition() != "North"){
		utility_exit_with_message("Antibody Design Protocol requires AHO_scheme and North definitions");
	}
	
	set_use_cluster_constraints(true);
	scorefxn_ = core::scoring::getScoreFunction();
	set_cdr_range(CDRNameEnum_start, CDRNameEnum_total, true);
	set_use_conservative_design_range(CDRNameEnum_start, CDRNameEnum_total, false);
	read_command_line_options();
	read_instructions(instruction_path_);//Any settings found in this file override cmdline defaults
	
}

AntibodyCDRDesigner::AntibodyCDRDesigner(AntibodyInfoOP ab_info, std::string instruction_path){
	ab_info_ = ab_info;
	if (ab_info_->get_current_AntibodyNumberingScheme()!="AHO_Scheme" && ab_info_->get_current_CDRDefinition() != "North"){
		utility_exit_with_message("Antibody Design Protocol requires AHO_scheme and North definitions");
	}
	
	set_use_cluster_constraints(true);
	scorefxn_ = core::scoring::getScoreFunction();
	set_cdr_range(CDRNameEnum_start, CDRNameEnum_total, true);
	set_use_conservative_design_range(CDRNameEnum_start, CDRNameEnum_total, false);
	read_command_line_options();
	instruction_path_=instruction_path;
	read_instructions(instruction_path_);//Any settings found in this file override cmdline defaults
}

AntibodyCDRDesigner::~AntibodyCDRDesigner(){}

void
AntibodyCDRDesigner::read_command_line_options(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	
	set_neighbor_detection_dis(option [OptionKeys::antibody::design::neighbor_dis]());
	set_use_conservative_design(h3, option [OptionKeys::antibody::design::conservative_h3_design]());
	set_use_turn_conservation(option [OptionKeys::antibody::design::turn_conservation]());
	set_probability_data_cutoff(option [OptionKeys::antibody::design::stats_cutoff]());
	set_basic_design(option [OptionKeys::antibody::design::benchmark_basic_design]());
	set_design_method(design_type_from_string(option [OptionKeys::antibody::design::design_method]()));
	set_rounds(option [OptionKeys::antibody::design::design_rounds]());
	instruction_path_ = basic::options::option [basic::options::OptionKeys::antibody::design::instructions]();
	zero_prob_weight_ = basic::options::option [basic::options::OptionKeys::antibody::design::sample_zero_probs_at]();
	
	
}

void
AntibodyCDRDesigner::read_instructions(std::string instruction_path){
	DesignInstructionsParser parser = DesignInstructionsParser(ab_info_, instruction_path);
	parser.read_cdr_design_instructions(instructions_);
}

std::string
AntibodyCDRDesigner::get_name() const {
	return "AntibodyCDRDesigner";
}

void
AntibodyCDRDesigner::set_cdr(const CDRNameEnum cdr, bool const setting){
	instructions_[cdr].design = setting;
}

void
AntibodyCDRDesigner::set_cdr_only(const CDRNameEnum cdr, bool const setting){
	if (setting==true){
		set_cdr_range(CDRNameEnum_start, CDRNameEnum_total, false);
		set_cdr(cdr, true);
	}
	else{
		set_cdr_range(CDRNameEnum_start, CDRNameEnum_total, true);
		set_cdr(cdr, false);
	}	
}

void
AntibodyCDRDesigner::set_cdr_range(const CDRNameEnum cdr_start, const CDRNameEnum cdr_end, bool const setting){
	for (core::SSize i=cdr_start; i<=cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		set_cdr(cdr, setting);
	}
}

void
AntibodyCDRDesigner::set_use_conservative_design(CDRNameEnum const cdr, bool const setting){
	instructions_[cdr].conservative_design = setting;
}

void
AntibodyCDRDesigner::set_use_conservative_design_range(const CDRNameEnum cdr_start, const CDRNameEnum cdr_end, const bool setting){
	for (core::SSize i=cdr_start; i<=cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		set_use_conservative_design(cdr, setting);
	}
}

void
AntibodyCDRDesigner::set_use_turn_conservation(const bool setting){
	turn_conservation_ = setting;
}

void
AntibodyCDRDesigner::set_basic_design(const bool setting){
	basic_design_ = setting;
}

void
AntibodyCDRDesigner::set_neighbor_detection_dis(core::Real const neighbor_distance){
	neighbor_dis_ = neighbor_distance;
}

void
AntibodyCDRDesigner::set_scorefxn(core::scoring::ScoreFunctionCOP scorefxn) {
	scorefxn_ = scorefxn->clone();
}

void
AntibodyCDRDesigner::set_rounds(const core::Size rounds){
	rounds_ = rounds;
}

void
AntibodyCDRDesigner::set_probability_data_cutoff(core::Size const cutoff){
	prob_cutoff_ = cutoff;
}

void
AntibodyCDRDesigner::set_zero_prob_weight_at(const core::Real weight){
	zero_prob_weight_ = weight;
}

void
AntibodyCDRDesigner::no_design_proline(const bool setting) {
	no_design_proline_ = setting;
}

void
AntibodyCDRDesigner::set_design_method(DesignTypeEnum const design_method){
	design_method_ = design_method;
}

void
AntibodyCDRDesigner::set_use_cluster_constraints(const bool setting) {
	use_cluster_constraints_ = setting;
}

std::map< core::Size, std::map< core::chemical::AA, core::Real > >
AntibodyCDRDesigner::setup_probability_data(core::pose::Pose& pose){
	AntibodyDatabaseManager manager = AntibodyDatabaseManager();
	std::map< core::Size, std::map< core::chemical::AA, core::Real > > prob_set;
	vector1<CDRNameEnum> no_data_cdrs= manager.load_cdr_design_data(ab_info_, pose, prob_set, prob_cutoff_, instructions_);
	
	for (core::Size i = 1; i<=no_data_cdrs.size(); ++i){
		instructions_[no_data_cdrs[i]].conservative_design = true;
	}
	
	return prob_set;
}

void
AntibodyCDRDesigner::remove_conservative_design_residues_from_prob_set(vector1<core::Size> const & positions, std::map< core::Size, std::map< core::chemical::AA, core::Real > > & prob_set){
	for (core::Size i = 1; i<= positions.size(); ++i){
		std::map< core::Size, std::map<core::chemical::AA, core::Real > >::iterator it = prob_set.find(positions[i]);
		if ( it != prob_set.end()){
			prob_set.erase(it);
			TR << "Removing "<<i << "  from probabilistic design" <<std::endl;
		}
	}
}

vector1<core::Size>
AntibodyCDRDesigner::get_conservative_design_residues(core::pose::Pose& pose){
	
	vector1<core::Size> conservative_positions;
	
	for (core::SSize i = CDRNameEnum_start; i <= CDRNameEnum_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		if (! instructions_[cdr].conservative_design) continue;
		
		core::Size start = ab_info_->get_CDR_start(cdr, pose);
		core::Size end = ab_info_->get_CDR_end(cdr, pose);
		
		for (core::Size res = start; res <= end; ++res){
			conservative_positions.push_back(res);
			//TR << "Treating "<< res << " as conservative " << std::endl;
		}
	}
	return conservative_positions;
}

void
AntibodyCDRDesigner::disable_design_cdrs(core::pack::task::TaskFactoryOP tf, core::pose::Pose & pose ) {
	for (core::Size i = 1; i <= CDRNameEnum_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		if (! instructions_[cdr].design){
			disable_design_cdr(cdr, tf, pose);
		}
	}
}

void
AntibodyCDRDesigner::disable_design_cdr(CDRNameEnum cdr, core::pack::task::TaskFactoryOP tf, core::pose::Pose & pose) {
	
	//One restrict op per CDR.  That way we can pop them off the TF  individually if we need to.
	core::pack::task::operation::RestrictResidueToRepackingOP restrict = new core::pack::task::operation::RestrictResidueToRepacking();
	core::Size start = ab_info_->get_CDR_start(cdr, pose);
	core::Size end = ab_info_->get_CDR_end(cdr, pose);
	for (core::Size i = start; i <= end; ++i){
		restrict->include_residue(i);
	}
	tf->push_back(restrict);
}

core::pack::task::TaskFactoryOP
AntibodyCDRDesigner::setup_task_factory(core::pose::Pose & pose){
	

	protocols::loops::LoopsOP cdr_loops = new protocols::loops::Loops();
	for (core::Size i = 1; i <= CDRNameEnum_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		if (instructions_[cdr].design){
			protocols::loops::Loop cdr_loop = ab_info_->get_CDR_loop(cdr, pose);
			cdr_loops->add_loop(cdr_loop);
		}
	}
	
	core::pack::task::TaskFactoryOP tf = new TaskFactory();
	
	//Setup Basic TaskOP
	tf->push_back(new InitializeFromCommandline());
	//tf->push_back(new RestrictToRepacking());
	
	//Setup Loops TaskOp + Turn on design for CDRs
	RestrictToLoopsAndNeighborsOP loop_task = new RestrictToLoopsAndNeighbors();
	loop_task->set_loops(cdr_loops);
	loop_task->set_design_loop(true);
	loop_task->set_include_neighbors(true);
	loop_task->set_cutoff_distance(neighbor_dis_);
	tf->push_back(loop_task);
	disable_design_cdrs(tf, pose);
	tf->push_back(new operation::NoRepackDisulfides());
	
	//Optionally disable Proline design
	
	//For benchmarking
	if (basic_design_){
		return tf;
	}
	
	//Setup Prob TaskOp.
	TR << "Adding ResidueProbDesignOp " << std::endl;
	ResidueProbDesignOperationOP prob_task = new ResidueProbDesignOperation();
	std::map< core::Size, std::map< core::chemical::AA, core::Real > > prob_set = setup_probability_data(pose);
	vector1<core::Size> conservative_positions = get_conservative_design_residues(pose);
	
	remove_conservative_design_residues_from_prob_set(conservative_positions, prob_set);
	
	prob_task->set_aa_probability_set(prob_set);
	prob_task->set_keep_task_allowed_aas(false);
	prob_task->set_include_native_restype(true);
	prob_task->set_sample_zero_probs_at(zero_prob_weight_); //Control by cmd line flag. We do want to have some variability that isn't known.;
	tf->push_back(prob_task);
	
	//Use conservative mutations for non-cluster positions + Optionally H3.
	
	if (!conservative_positions.empty()){
		TR << "Adding ConservativeDesignOp "<<std::endl;
		ConservativeDesignOperationOP cons_task = new ConservativeDesignOperation();
		cons_task->limit_to_positions(conservative_positions);
		cons_task->include_native_aa(true);
		cons_task->add_to_allowed_aas(false);
		tf->push_back(cons_task);
	}
	
	return tf;
}

bool
AntibodyCDRDesigner::cdr_has_constraints(core::pose::Pose const & pose, CDRNameEnum const cdr, std::string const constraint_type){
	using namespace core::scoring::constraints;
	
	core::Size start_res = ab_info_->get_CDR_start(cdr, pose);
	core::Size end_res = ab_info_->get_CDR_end(cdr, pose);
	
	std::map< core::Size, bool> cst_found;
	
	//Initialize our map of whether the CDR residue has constraints.  
	for (core::Size i = start_res; i <= end_res; ++i){
		cst_found[i] = false;
	}
	utility::vector1< ConstraintCOP > csts = pose.constraint_set()->get_all_constraints();
	for (core::Size i = 1; i <= csts.size(); ++i){
		if (csts[i]->type() != constraint_type){ continue; }
		
		utility::vector1< core::Size > residues = csts[i]->residues();
		for (core::Size x = 1; x <= residues.size(); ++x){
			cst_found[residues[x]] = true;
		}
	}
	
	//Check that all residues have constraints of the particular type:
	for (core::Size i = start_res ; i <= end_res; ++i){
		if (! cst_found[i]){ return false; }
	}
	
	return true;
}

void
AntibodyCDRDesigner::setup_constraints(core::pose::Pose & pose){
	
	//We only add dihedral or coordinate constraints if they are not present for the pose at hand in each CDR.
	for (core::SSize i = 1; i<= CDRNameEnum_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		bool constraint_result = false;
		if (! cdr_has_constraints(pose, cdr, "DihedralConstraint")){
			constraint_result = protocols::antibody::add_harmonic_cluster_constraint(ab_info_, pose, ab_info_->get_CDR_cluster(cdr)->cluster());
		}
		if (! constraint_result && ! cdr_has_constraints(pose, cdr, "CoordinateConstraint")){
			core::Size start_res = ab_info_->get_CDR_start(cdr, pose);
			core::Size end_res = ab_info_->get_CDR_end(cdr, pose);
			//This needs to change to dihedral constraints - but what will the standard deviation be?
			
			TR << "Adding coordinate constraints for " << ab_info_->get_CDR_name(cdr) << std::endl;
			core::scoring::constraints::add_coordinate_constraints(pose, start_res, end_res, .5, false /* include_sc */);
		}
	}
}


void
AntibodyCDRDesigner::apply(core::pose::Pose& pose){

		using namespace protocols::antibody;
		using namespace protocols::antibody::design;
		using namespace core::pack::task;
		using namespace protocols::toolbox::task_operations;

		if (! protocols::antibody::clusters::check_if_pose_renumbered_for_clusters(pose)){
			utility_exit_with_message("PDB must be numbered correctly to identify North CDR clusters.  Please see Antibody Design documentation.");
		}
		
		ab_info_->setup_CDR_clusters(pose);
		core::pack::task::TaskFactoryOP tf = setup_task_factory(pose);
		

		
		protocols::simple_moves::PackRotamersMover packer = protocols::simple_moves::PackRotamersMover(scorefxn_);
		packer.task_factory(tf);
		
		//////////Setup//////////////////////////

		
		scorefxn_->set_weight(chainbreak, 100);
		
		if (use_cluster_constraints_){
			if (scorefxn_->get_weight(dihedral_constraint) == 0.0){
				scorefxn_->set_weight(dihedral_constraint, 1.0);
			}
			setup_constraints(pose);
		}
		
		
		scorefxn_->show(pose);
		
		//Setup the movemap.  Allow SC minimization with neighbors.
		core::kinematics::MoveMapOP movemap = new MoveMap();
		for (core::Size i = 1; i<=CDRNameEnum_total; ++i){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			if (instructions_[cdr].design){
				ab_info_->add_CDR_to_MoveMap(pose, movemap, cdr, false, true, neighbor_dis_);
			}
		}
		
		//Setup MonteCarlo for >1 round
		protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo(pose, *scorefxn_, 1.0);
		
		//Design methods
		if (design_method_ == fixbb){
			for (core::Size i = 1; i <= rounds_; ++i){
				packer.apply(pose);
				TR << (*scorefxn_)(pose) <<std::endl;
				mc->boltzmann(pose);
			}
		}
		else if( design_method_ == flxbb) {
			protocols::flxbb::FlxbbDesignOP flx = new protocols::flxbb::FlxbbDesign();
			protocols::flxbb::DesignTaskOP des = new protocols::flxbb::DesignTask_Normal();
			flx->set_movemap(movemap);
			core::pack::task::TaskFactory::const_iterator it;
			for (it = tf->begin(); it != tf->end(); ++it){
				des->add_task_operation(*it);
			}
			
			flx->add_design_task(des);
			flx->set_scorefxn_design(scorefxn_);
			flx->set_scorefxn_relax(scorefxn_);
			
			for (core::Size i =1; i <= rounds_; ++i){
				flx->apply(pose);
				TR << (*scorefxn_)(pose) <<std::endl;
				mc->boltzmann(pose);
			}
			
			
		}
		else if ( design_method_ == relaxed_design) {

			
			protocols::relax::FastRelaxOP rel = new protocols::relax::FastRelax(scorefxn_);
			rel->set_movemap(movemap);
			
			//Optionally minimize bond length + angles - Test using cmd-line first!
			rel->set_task_factory(tf);
			for (core::Size i = 1; i<=rounds_; ++i){
				rel->apply(pose);
				TR << (*scorefxn_)(pose) <<std::endl;
				mc->boltzmann(pose);
			}
			
		}
		else if ( design_method_ == docked_design){
			TR << "Docked design not currently implemented. " << std::endl;
		}
		else{
			utility_exit_with_message("Antibody Design method not recognized. ");
		}
		
		
		mc->recover_low(pose);
		TR << "AntibodyCDRDesigner complete." << std::endl;
}





} //design
} //antibody
} //protocols
