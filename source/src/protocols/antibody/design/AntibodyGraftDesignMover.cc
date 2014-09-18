// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/design/AntibodyGraftDesignMover.cc
/// @brief Class that initially designs antibodies through grafting using an AntibodyDatabase + North_AHO numbering scheme
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Project Includes
#include <protocols/antibody/design/AntibodyGraftDesignMover.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>
#include <protocols/antibody/clusters/util.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/util.hh>

#include <protocols/antibody/design/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/DesignInstructionsParser.hh>
#include <protocols/antibody/design/AntibodyDesignModeler.hh>
#include <protocols/antibody/design/util.hh>

// Core Includes
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/util.hh>

#include <core/pose/PDBInfo.hh>

// Protocol Includes
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/grafting/CCDEndsGraftMover.hh>
#include <protocols/grafting/util.hh>
#include <protocols/forge/methods/chainbreak_eval.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/util.hh>

// Numeric Includes
#include <numeric/random/random.hh>

// Options
#include <basic/options/option.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

//Utility
#include <math.h>
#include <basic/Tracer.hh>
#include <utility/PyAssert.hh>
#include <map>
#include <boost/algorithm/string.hpp>


static thread_local basic::Tracer TR( "antibody.design.AntibodyGraftDesinger" );

namespace protocols{
namespace antibody{
namespace design{
	using namespace protocols::antibody;
	using namespace protocols::grafting;
	using namespace protocols::antibody::clusters;
	using namespace core;
	using core::Size;

AntibodyGraftDesignMover::AntibodyGraftDesignMover(AntibodyInfoOP ab_info):
	graft_mover_(NULL),
	scorefxn_(NULL)
{
	overhang_ = 3;
	ab_info_=ab_info;
	max_linear_chainbreak_ = .6;
	modeler_ = new AntibodyDesignModeler(ab_info_);
	if (ab_info_->get_current_AntibodyNumberingScheme()!="AHO_Scheme" && ab_info_->get_current_CDRDefinition() != "North"){
		utility_exit_with_message("Antibody Design Protocol requires AHO_scheme and North definitions");
	}
	read_command_line_options();
	set_defaults();
}

AntibodyGraftDesignMover::AntibodyGraftDesignMover(AntibodyInfoOP ab_info, std::string instruction_path) :
	graft_mover_(NULL),
	scorefxn_(NULL)
{
	overhang_ = 3;
	ab_info_=ab_info;
	max_linear_chainbreak_ = .35;
	modeler_ = new AntibodyDesignModeler(ab_info_);
	if (ab_info_->get_current_AntibodyNumberingScheme()!="AHO_Scheme" && ab_info_->get_current_CDRDefinition() != "North"){
		utility_exit_with_message("Antibody Design Protocol requires AHO_scheme and North definitions");
	}
	read_command_line_options();
	instruction_path_ = instruction_path;
	set_defaults();
}

AntibodyGraftDesignMover::~AntibodyGraftDesignMover(){}

void
AntibodyGraftDesignMover::set_defaults(){
	///Conservative defaults.  Defaults here are also read and set from a database file.  To allow run-time manipulations and testing.

	set_cdr_range(CDRNameEnum_start, CDRNameEnum_total, true);
	set_mintype_range(CDRNameEnum_start, CDRNameEnum_total, minimize);
	set_min_rb_range(CDRNameEnum_start, CDRNameEnum_total, false);
	set_min_neighbor_sc_range(CDRNameEnum_start, CDRNameEnum_total, true);
	set_cdr_range_stay_in_native_cluster(CDRNameEnum_start, CDRNameEnum_total, false);
	set_cdr_range_stay_in_type(CDRNameEnum_start, CDRNameEnum_total, 1, true);
	set_cdr_cluster_centers_only_range(CDRNameEnum_start, CDRNameEnum_total, false);
	set_cdr_min_length_range(CDRNameEnum_start, CDRNameEnum_total, 1);
	set_cdr_max_length_range(CDRNameEnum_start, CDRNameEnum_total, 50); //Something absurd.
	read_instructions(instruction_path_);//Any settings found in this file override these defaults.
}

void
AntibodyGraftDesignMover::read_command_line_options(){
	set_graft_rounds(basic::options::option [basic::options::OptionKeys::antibody::design::graft_rounds]());
	instruction_path_ = basic::options::option [basic::options::OptionKeys::antibody::design::instructions]();
	set_keep_top_designs(basic::options::option [basic::options::OptionKeys::antibody::design::top_graft_designs]());
	set_dock_post_graft(basic::options::option [basic::options::OptionKeys::antibody::design::dock_post_graft]());
	set_pack_post_graft(basic::options::option [basic::options::OptionKeys::antibody::design::pack_post_graft]());
	set_rb_min_post_graft(basic::options::option [basic::options::OptionKeys::antibody::design::rb_min_post_graft]());
	set_dock_rounds(basic::options::option [basic::options::OptionKeys::antibody::design::dock_rounds]());
	initial_perturb_ = basic::options::option [basic::options::OptionKeys::antibody::design::initial_perturb] ();
	use_deterministic_algorithm_ = basic::options::option [basic::options::OptionKeys::antibody::design::use_deterministic]();
	//benchmark_ = basic::options::option [basic::options::OptionKeys::antibody::design::benchmark_graft_designer]();
	extend_native_cdrs_ = basic::options::option [basic::options::OptionKeys::antibody::design::extend_native_cdrs]();
}

void
AntibodyGraftDesignMover::set_scorefunction(ScoreFunctionOP scorefxn){
	scorefxn_=scorefxn->clone();
	if (scorefxn_->get_weight(core::scoring::dihedral_constraint)==0){
		scorefxn_->set_weight(core::scoring::dihedral_constraint, 1.0);
	}
}

void
AntibodyGraftDesignMover::set_keep_top_designs(core::Size top_designs){
	num_top_designs_ = top_designs;
}

void
AntibodyGraftDesignMover::set_dock_post_graft(bool dock_post_graft){
	dock_post_graft_ = dock_post_graft;
}

void
AntibodyGraftDesignMover::set_dock_rounds(core::Size dock_rounds){
	dock_rounds_ = dock_rounds;
}

void
AntibodyGraftDesignMover::set_pack_post_graft(bool pack_post_graft){
	pack_post_graft_ = pack_post_graft;
}

void
AntibodyGraftDesignMover::set_rb_min_post_graft(bool rb_min_post_graft){
	rb_min_post_graft_ = rb_min_post_graft;
}
void
AntibodyGraftDesignMover::set_graft_rounds(core::Size graft_rounds){
	graft_rounds_=graft_rounds;
}

void
AntibodyGraftDesignMover::set_cluster(const CDRNameEnum cdr, const CDRClusterEnum cluster, const bool setting){
	if (setting==true){
		cdr_instructions_[cdr].leave_out_clusters.push_back(cluster);
	}
	else{
		cdr_instructions_[cdr].include_only_clusters.push_back(cluster);
	}
}

void
AntibodyGraftDesignMover::set_cdr_set(CDRSet& cdr_set, CDRClusterMap & cdr_cluster_map, core::Size overhang){
	cdr_set_ = cdr_set;
	cdr_cluster_map_ = cdr_cluster_map;
	overhang_ = overhang;
}

void
AntibodyGraftDesignMover::set_cluster_range(
		const CDRNameEnum cdr,
		const CDRClusterEnum cluster_start,
		const CDRClusterEnum cluster_end,
		const bool setting)
{
	for (core::SSize i=cluster_start; i<=cluster_end; ++i){
		CDRClusterEnum cluster = static_cast<CDRClusterEnum>(i);
		if (setting==true){

			cdr_instructions_[cdr].leave_out_clusters.push_back(cluster);
		}
		else{
			cdr_instructions_[cdr].include_only_clusters.push_back(cluster);
		}
	}
}

void
AntibodyGraftDesignMover::set_mintype(const CDRNameEnum cdr_name, MinTypeEnum mintype){
	cdr_instructions_[cdr_name].mintype = mintype;
}

void
AntibodyGraftDesignMover::set_mintype_range(const CDRNameEnum cdr_start, const CDRNameEnum cdr_end, MinTypeEnum mintype){
	for (core::SSize i=cdr_start; i<=cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		set_mintype(cdr, mintype);
	}
}

void
AntibodyGraftDesignMover::set_min_rb(const CDRNameEnum cdr_name, const bool setting){
	cdr_instructions_[cdr_name].min_rb = setting;
}

void
AntibodyGraftDesignMover::set_min_rb_range(const CDRNameEnum cdr_start, const CDRNameEnum cdr_end, bool const setting){
	for (core::SSize i = cdr_start; i <= cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		set_min_rb(cdr, setting);
	}
}
void
AntibodyGraftDesignMover::set_min_neighbor_sc(const CDRNameEnum cdr_name, bool const setting){
	cdr_instructions_[cdr_name].min_neighbor_sc = setting;
}

void
AntibodyGraftDesignMover::set_min_neighbor_sc_range(const CDRNameEnum cdr_start, const CDRNameEnum cdr_end, bool const setting){
	for (core::SSize i = cdr_start; i <= cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		set_min_neighbor_sc(cdr, setting);
	}
}

void
AntibodyGraftDesignMover::set_cdr(const CDRNameEnum cdr, const bool setting){
	//cdrs_to_use_[cdr]=setting;
	cdr_instructions_[cdr].graft= setting;
}

void
AntibodyGraftDesignMover::set_cdr_range(const CDRNameEnum cdr_start, const CDRNameEnum cdr_end, const bool setting){
	for (core::SSize i=cdr_start; i<=cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		cdr_instructions_[cdr].graft = setting;
	}
}

void
AntibodyGraftDesignMover::set_cdr_stay_in_native_cluster(const CDRNameEnum cdr, const bool setting){
	cdr_instructions_[cdr].stay_native_cluster = setting;
}

void
AntibodyGraftDesignMover::set_cdr_range_stay_in_native_cluster(

				const CDRNameEnum cdr_start,
				const CDRNameEnum cdr_end,
				const bool setting)
{

	for (core::SSize i=cdr_start; i<=cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		cdr_instructions_[cdr].stay_native_cluster = setting;
	}
}

void
AntibodyGraftDesignMover::set_cdr_stay_in_type(const CDRNameEnum cdr, const core::Size type, const bool setting){
	cdr_instructions_[cdr].cluster_types[type] = setting;
}

void
AntibodyGraftDesignMover::set_cdr_range_stay_in_type(

				const CDRNameEnum cdr_start,
				const CDRNameEnum cdr_end,
				const core::Size type,
				const bool setting)
{
	for (core::SSize i=cdr_start; i<=cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		cdr_instructions_[cdr].cluster_types[type]=setting;
	}
}

void
AntibodyGraftDesignMover::set_cdr_min_length(CDRNameEnum const cdr, core::Size length){
	cdr_instructions_[cdr].min_length = length;
}

void
AntibodyGraftDesignMover::set_cdr_min_length_range(const CDRNameEnum cdr_start, CDRNameEnum cdr_end, core::Size length){
	for (core::SSize i=cdr_start; i<=cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		set_cdr_min_length(cdr, length);
	}
}

void
AntibodyGraftDesignMover::set_cdr_max_length(CDRNameEnum const cdr, core::Size length){
	cdr_instructions_[cdr].max_length = length;
}

void
AntibodyGraftDesignMover::set_cdr_max_length_range(const CDRNameEnum cdr_start,const CDRNameEnum cdr_end, core::Size length){
	for (core::SSize i=cdr_start; i<=cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		set_cdr_max_length(cdr, length);
	}
}

void
AntibodyGraftDesignMover::set_cdr_cluster_centers_only(const CDRNameEnum cdr, bool setting){
	cdr_instructions_[cdr].cluster_centers_only = setting;
}

void
AntibodyGraftDesignMover::set_cdr_cluster_centers_only_range(const CDRNameEnum cdr_start, CDRNameEnum cdr_end, bool setting){
	for (core::SSize i=cdr_start; i<=cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		set_cdr_cluster_centers_only(cdr, setting);
	}
}

void
AntibodyGraftDesignMover::setup_native_clusters(pose::Pose & pose){

	ab_info_->setup_CDR_clusters(pose);

	for (Size i=1; i<=CDRNameEnum_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		native_clusters_[cdr] = ab_info_->get_CDR_cluster(cdr);
	}
}

void
AntibodyGraftDesignMover::initialize_cdr_set(){
	AntibodyDatabaseManagerOP manager = new AntibodyDatabaseManager();
	std::pair< CDRSet,  CDRClusterMap > result_pair = manager->load_cdrs_for_grafting(ab_info_, cdr_instructions_, pdbmap_);
	cdr_set_ = result_pair.first; cdr_cluster_map_ = result_pair.second;

}

void
AntibodyGraftDesignMover::read_instructions(std::string instruction_path){
	DesignInstructionsParser parser = DesignInstructionsParser(ab_info_, instruction_path);
	parser.read_cdr_graft_instructions(cdr_instructions_);
}

void
AntibodyGraftDesignMover::fix_pdb_info(pose::Pose & pose, CDRNameEnum cdr, CDRClusterEnum cluster, core::Size original_start, core::Size original_pdb_end) {

	//using core::pose::PDBInfo;

	//PDBInfo should not be obsolete when calling this function.
	//PDBInfo & pdbinfo = *pose.pdb_info();
	core::Size cdr_length = ab_info_->get_cluster_length(cluster);
	core::Real result = (core::Real)cdr_length/2;
	core::Size pdb_num_start = pose.pdb_info()->number(original_start);
	core::Size pdb_num_end = original_pdb_end;
	core::Size original_end = original_start+cdr_length+1;

	//Starting residues
	for (core::Size i = 1; i<=ceil(result); ++i){
		core::Size pdb_num = pdb_num_start+i;
		core::Size pose_num = original_start+i;
		char chain = ab_info_->get_CDR_chain(cdr);
		TR.Debug << "Set PDBInfo: "<< pose_num<<": "<< pdb_num << " "<<  chain <<std::endl;
		pose.pdb_info()->set_resinfo(pose_num, chain , pdb_num);
	}

	//Ending residues
	for (core::Size i = 1; i<=floor(result); ++i){
		core::Size pdb_num = pdb_num_end-i;
		core::Size pose_num = original_end-i;
		char chain = ab_info_->get_CDR_chain(cdr);
		TR.Debug << "Set PDBInfo: "<< pose_num<<": "<< pdb_num << " "<<  chain <<std::endl;
		pose.pdb_info()->set_resinfo(pose_num, chain , pdb_num);
	}
}

void
AntibodyGraftDesignMover::set_default_graft_settings(){
	graft_mover_->set_cycles(20);
	graft_mover_->set_scaffold_flexibility(2, 2);
	graft_mover_->set_insert_flexibility(0, 0);
	graft_mover_->final_repack(true);
	graft_mover_->stop_at_closure(true);
	graft_mover_->neighbor_dis(modeler_->get_neighbor_detection_dis());
	//graft_mover_->set_skip_sampling(true);
}

void
AntibodyGraftDesignMover::extend_native_cdrs(pose::Pose & pose, vector1<CDRNameEnum>& cdrs_to_design) {

	for (core::Size i = 1; i <= cdrs_to_design.size(); ++i){
		TR << "Randomizing native CDR structure: "<< ab_info_->get_CDR_name(cdrs_to_design[i]) << std::endl;
		modeler_->extend_CDR(pose, cdrs_to_design[i]);
	}
	//pose.dump_pdb("extend_test.pdb");
}

bool
AntibodyGraftDesignMover::graft_cdr(pose::Pose& pose, CDRNameEnum cdr, core::Size index){

	//Get CDR start/end in PDBInfo is now dynamic.  As long as pdb_info is not obsolete
	//This should go into a log, we should figure out how to have this info written to the resulting file.
	TR << "Grafting CDR from cluster " << ab_info_->get_cluster_name(cdr_cluster_map_[cdr][index]) << " fragment "<< pdbmap_[cdr][index] << std::endl;

	core::Size start = ab_info_->get_CDR_start(cdr, pose)-1;
	core::Size end = ab_info_->get_CDR_end(cdr, pose)+1;
	//core::Size original_pdb_end = pose.pdb_info()->number(end);

	graft_mover_->set_insert_region(start, end);
	//(cdr_set_[cdr][index])->dump_pdb("piece.pdb");

	if (cdr==h3){
		//graft_mover_->set_scaffold_flexibility(3, 3);
		graft_mover_->set_cycles(40);
	}

	//If the superposition is really wacky for some reason, we can get chainbreaks.  To overcome this,:
	//Next - need to min on a jump or multiple jumps as this adaptation does not always work - but it will take too many cycles to use normal graft algorithm.
	core::Size max_attempts = 1;
	pose::Pose temp_pose = pose;
	for (core::Size i = 1; i <= max_attempts; ++i){
		pose::Pose piece = pose::Pose(*(cdr_set_[cdr][index]));
		graft_mover_->set_piece(piece, overhang_, overhang_);
		//graft_mover_->superimpose_overhangs_heavy(temp_pose, ca_only , true); //Superimposes first heavy atoms of overhang
		graft_mover_->apply(temp_pose);


		//Is there a chainbreak?
		bool cb = false;
		core::Real cb_score;
		TR << "Checking chainbreak" << std::endl;
		for (core::Size pos = start - graft_mover_->get_nterm_scaffold_flexibility(); pos <= end + graft_mover_->get_cterm_scaffold_flexibility(); ++pos){
			cb_score = protocols::forge::methods::linear_chainbreak(temp_pose, pos);
			if (cb_score >= max_linear_chainbreak_){
				TR << "Chainbreak found" << std::endl;
				cb = true;
				break;
			}
		}

		//Increase sampling to do a better graft
		if (cb && i <  max_attempts){
			TR << "Linear chainbreak found from graft. Increasing flexibility." << std::endl;
			TR << "Value: " << cb_score << std::endl;
			//graft_mover_->set_use_double_loop_double_CCD_arms(false);
			//graft_mover_->set_use_double_loop_quad_CCD_arms(true);
			graft_mover_->set_scaffold_flexibility(graft_mover_->get_nterm_scaffold_flexibility()+1, graft_mover_->get_cterm_scaffold_flexibility()+1);
			graft_mover_->set_insert_flexibility(1, 1);
			//graft_mover_->set_cycles(20);
			//ca_only = true;
		}
		else if (cb && i == max_attempts){
			//Still a chainbreak
			set_default_graft_settings();
			TR << "Could not close graft for CDR after adaptation...  Skipping." << std::endl;
			TR << "Value: " << cb_score << std::endl;
			return false;
		}
		else{
			//Success
			pose = temp_pose;
			set_default_graft_settings();
			break;
		}
	} // End graft adaptation

	pose.pdb_info()->copy(*((cdr_set_[cdr][index])->pdb_info()), 1 + overhang_, (cdr_set_[cdr][index])->total_residue()-overhang_, start+1);
	pose.pdb_info()->obsolete(false);

	TR << "PDBInfo Set" << std::endl;
	modeler_->set_cdr_only(cdr, true);


	//Remove old CDR constraints.  (No need to make local copy due to MC holding the lowest pose - if we go back, constraints will be there.)
	pose.remove_constraints(constraint_map_[cdr], true);
	core::Size start_res = ab_info_->get_CDR_start(cdr, pose);
	core::Size end_res = ab_info_->get_CDR_end(cdr, pose);
	core::scoring::constraints::remove_constraints_of_type(pose, "CoordinateConstraint", start_res, end_res);

	//Add new post-graft CDR constraints or coordinate constraints if unknown cluster.
	bool constraint_result = protocols::antibody::add_harmonic_cluster_constraint(ab_info_, pose, cdr_cluster_map_[cdr][index], constraint_map_[cdr]);
	bool use_coordinate_constraints = !constraint_result;
	if (use_coordinate_constraints){

		//This needs to change to dihedral constraints - but what will the standard deviation be?
		TR << "Adding coordinate constraints for "<< ab_info_->get_CDR_name(cdr) << std::endl;
		core::scoring::constraints::add_coordinate_constraints(pose, start_res, end_res, .5 /* just a little give */, false /* include_sc */);
	}


	///Relax CDR if in instructions + the constraint was successfully added.
	vector1< char > antigen_chains = ab_info_->get_antigen_chains();
	std::string antigen(antigen_chains.begin(), antigen_chains.end());
	std::string dock_chains = modeler_->get_ab_dock_chains()+"_"+antigen;

	TR << "Beginning post-graft min" << std::endl;
	if (initial_perturb_){
		core::kinematics::FoldTree ft = pose.fold_tree();
		vector1< int > movable_jumps(1, 1);
		protocols::docking::setup_foldtree(pose, dock_chains, movable_jumps);
		protocols::docking::DockingInitialPerturbationOP perturber = new protocols::docking::DockingInitialPerturbation(1, true /* slide */);

		perturber->apply(pose);
		pose.fold_tree(ft);
	}

	if (pack_post_graft_ && cdr_instructions_[cdr].mintype != no_min){
		modeler_->repack_CDRs_and_neighbors(pose);
	}

	if (dock_post_graft_){
		for (core::Size i = 1; i<= dock_rounds_; i++){

			//Control this from cmd line
			//Perturb if command-line options are set: At the very least, we slide into contact with the antibody.


			modeler_->dock_LH_A_low_res(pose, false /*repack_interface*/ ); // This should change once the minimization does some neighbor detection.
			modeler_->dock_LH_A_high_res(pose, 3 /*first cycles*/, 10 /*second_cycles*/); //Normal DockMCM is 4/45.  This should mainly just be quick to fix overlap from low_res dock.
		}
	}



	bool rb_min = cdr_instructions_[cdr].min_rb;

	//I would have liked to control this via a RosettaScript, read in during the run.  But, members in Sarels lab say this is not possible yet.
	if (cdr_instructions_[cdr].mintype == relax){

		if (cdr_instructions_[cdr].min_neighbor_sc){
			modeler_->relax_cdrs_and_neighbor_sc(pose, false, rb_min, dock_chains);
		}
		else{
			modeler_->relax_cdrs(pose, false /*centroid_mode*/ , false, rb_min, dock_chains);//Change this if need be.
		}
	}
	else if (cdr_instructions_[cdr].mintype == centroid_relax){
		modeler_->relax_cdrs(pose, true /*centroid_mode*/ , false, rb_min, dock_chains);//Change this if need be.
	}
	//Side chains are already repacked by default after graft mover
	else if (cdr_instructions_[cdr].mintype == minimize){
		//graft_mover_->repack_connection_and_residues_in_movemap_and_piece(pose, scorefxn_);
		if (cdr_instructions_[cdr].min_neighbor_sc){
			modeler_->minimize_cdrs_and_neighbor_sc(pose, rb_min, dock_chains);
		}
		else {
			modeler_->minimize_cdrs(pose, rb_min, dock_chains);
		}
	}
	else {
		//graft_mover_->repack_connection_and_residues_in_movemap(pose, scorefxn_);
	}

	if (rb_min_post_graft_){
		modeler_->minimize_interface(pose, dock_chains, true /* min_interface */);
	}
	TR << "Graft complete" << std::endl;
	return true;
}

void
AntibodyGraftDesignMover::check_for_top_designs(pose::Pose & pose){

	//Can be refactored to use utility::TopScoreSelector
	//From mc algorithm, you can have multiple poses that are equivalent...
	core::Real score = (*scorefxn_)(pose);


	vector1<core::Real>::iterator score_it = top_scores_.begin();
	vector1<pose::PoseOP>::iterator pose_it = top_designs_.begin();

	if (top_scores_.size()==0){
		top_scores_.push_back(score);
		top_designs_.push_back(new Pose());
		*(top_designs_[top_designs_.size()]) = pose;
	}
	else{
		bool inserted = false;
		for (core::Size i = 1; i<=top_scores_.size(); ++i){
			if (score <= top_scores_[i]){
				top_scores_.insert(score_it+i-1, score);
				top_designs_.insert(pose_it+i-1, new Pose());
				*(top_designs_[i]) = pose;
				inserted = true;
				break;
			}
		}
		if (! inserted && top_scores_.size() < num_top_designs_){
			top_scores_.push_back(score);
			top_designs_.push_back(new Pose());
			*(top_designs_[top_designs_.size()]) = pose;
		}
		else if ( inserted && top_scores_.size() > num_top_designs_){
			top_scores_.pop_back();
			top_designs_.pop_back();
		}
	}

	//mc_->eval_lowest_score_pose(pose, false, true);

}

///@brief Gets a list of vectors whose indexes correspond to CDRNameEnum, and whose values correspond to the cdr_set index.  If the value is 0, it means no cdr in set.
vector1< vector1 < Size > >
AntibodyGraftDesignMover::get_cdr_set_index_list(){
	vector1< vector1< Size > > index_list;
	vector1<core::Size> cdr_set_totals(6, 0);

	for (core::SSize i=CDRNameEnum_start; i<=CDRNameEnum_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		cdr_set_totals[i] = cdr_set_[cdr].size();
	}

	vector1< Size > dummy_index(6, 0);
	get_all_graft_permutations(cdr_set_totals, index_list, dummy_index, 1);
	return index_list;

}

void
AntibodyGraftDesignMover::run_basic_mc_algorithm(pose::Pose& pose, vector1<CDRNameEnum>& cdrs_to_design){
	TR << "Running basic monte carlo algorithm " << std::endl;

	top_scores_.push_back((*scorefxn_)(pose));
	top_designs_.push_back(new Pose());
	*(top_designs_[1]) = pose;
	mc_->set_last_accepted_pose(pose);

	for (core::Size i = 1; i <= graft_rounds_; ++i){
		TR << "Graft round: " << i <<std::endl;

		CDRNameEnum cdr_type = cdrs_to_design[numeric::random::rg().random_range(1, cdrs_to_design.size())];
		core::Size cdr_index = numeric::random::rg().random_range(1, cdr_set_[cdr_type].size());

		bool graft_successful = graft_cdr(pose, cdr_type, cdr_index);
		if ( graft_successful ){
			check_for_top_designs(pose);
			mc_->boltzmann(pose);
		}
		else{
			pose = mc_->last_accepted_pose();
		}
	}

	mc_->recover_low(pose);
	mc_->show_counters();
}

void
AntibodyGraftDesignMover::run_deterministic_graft_algorithm(pose::Pose & pose,vector1<CDRNameEnum>& cdrs_to_design){
	//Note:  Not feasible with >= 4 CDRs to try at each position. Should not be used with docking on.

	/*
	vector1< vector1< Size > > cdr_index_list = get_cdr_set_index_list();

	TR << "Running deterministic graft algorithm on  " << cdr_index_list.size() << " permutations. " << std::endl;
	for (core::Size i = 1; i <= cdr_index_list.size(); ++i){
		TR << "Graft round: " << i << std::endl;
		pose::Pose trial_pose = pose;

		//Graft each CDR in the index from the starting pose.  Same order every time.
		//All grafts should be closed for this to be OK.

		for (core::Size x = 1; x <= 6; ++x) {
			if (cdr_index_list[i][x] != 0){
				CDRNameEnum cdr = static_cast<CDRNameEnum>(x);
				bool graft_successful = graft_cdr(trial_pose, cdr, cdr_index_list[i][x]);
				if (! graft_successful){
					TR << "Graft unsuccessful.  Skipping combination..." << std::endl;
					break;
				}
			}
			//Only check for top designs after all combinations have been grafted.
			check_for_top_designs(trial_pose);
			mc_->eval_lowest_score_pose(trial_pose, false, true);
		}

	} //End grafting
	*/

	//Temporary fix to deterministically graft one CDR.  Needs to be fixed correctly soon.
	CDRNameEnum cdr = cdrs_to_design[1];

	for (core::Size i = 1; i <= cdr_set_[cdr].size(); ++i){
		TR << "Graft round: " << i << std::endl;
		pose::Pose trial_pose = pose;
		bool graft_successful = graft_cdr(trial_pose, cdr, i);
		if (! graft_successful){
			TR << "Graft unsuccessful.  Skipping combination..." << std::endl;
			continue;
		}
		//Only check for top designs after all combinations have been grafted.
		check_for_top_designs(trial_pose);
		mc_->eval_lowest_score_pose(trial_pose, false, true);
	}

	mc_->recover_low(pose);
	mc_->show_counters();
}

void
AntibodyGraftDesignMover::apply(pose::Pose & pose){

	if (! protocols::antibody::clusters::check_if_pose_renumbered_for_clusters(pose)){
		utility_exit_with_message("PDB must be numbered correctly to identify North CDR clusters.  Please see Antibody Design documentation.");
	}

	//Create Instances.  Make sure everything is ready to begin.
	if (!scorefxn_){
		scorefxn_ = core::scoring::get_score_function(true);
	}
	if (scorefxn_->get_weight(core::scoring::dihedral_constraint) == 0.0){
		scorefxn_->set_weight(core::scoring::dihedral_constraint, 1.0);
	}

	//Print setup options.  Will be moved to Show
	// -> Everything set
	// -> Total number of grafts
	// -> Total number of CDRs
	// -> Total number of CDR pieces per CDR

	//Will go in SHOW.
	TR << "///////////////////////////////////////////////////////////////////////////////////////////////////////////////" <<std::endl;
	for (core::Size i=1; i<=6; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		TR << "//////// " << ab_info_->get_CDR_name(cdr) << " ///////////////////////////////////////////////"<< std::endl;
		TR << "//// "<< std::endl;
		TR << "///  Graft? " << std::boolalpha << cdr_instructions_[cdr].graft << std::endl;
		//TR << "///  Relax? " << std::boolalpha << cdr_instructions_[cdr].relax << std::endl;
		TR<< "///  Stay Native? " << std::boolalpha << cdr_instructions_[cdr].stay_native_cluster << std::endl;
		TR<< "///  Cluster Centers only? " << std::boolalpha << cdr_instructions_[cdr].cluster_centers_only << std::endl;
		if (cdr!=h3){
			TR << "////" << std::endl;
			TR << "////// Types: " << std::endl;
			TR << "////" << std::endl;
			TR << "///  1 " << std::boolalpha << cdr_instructions_[cdr].cluster_types[1] << std::endl;
			TR << "///  2 " << std::boolalpha << cdr_instructions_[cdr].cluster_types[2] << std::endl;
			TR << "///  3 " << std::boolalpha << cdr_instructions_[cdr].cluster_types[3] << std::endl;
		}
		TR << "////" << std::endl;
		TR << "////// Lengths: " << std::endl;
		TR << "////" << std::endl;
		TR << "/// Min " << cdr_instructions_[cdr].min_length << std::endl;
		TR << "/// Max " << cdr_instructions_[cdr].max_length << std::endl;
		TR << "///" << std::endl;
	}
	TR << "///////////////////////////////////////////////////////////////////////////////////////////////////////////////" <<std::endl;
	setup_native_clusters(pose);
	if (cdr_set_.empty()){
		initialize_cdr_set();
	}

	//Reinitialize.  Need to figure out how to do this properly via JD.
	top_designs_.clear();
	top_scores_.clear();

	//List of CDRs to graft.
	vector1< CDRNameEnum > cdrs_to_design;

	for (core::Size i=1; i<=CDRNameEnum_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		if (cdr_instructions_[cdr].graft && cdr_set_[cdr].size() != 0){
			cdrs_to_design.push_back(cdr);
		}
	}
	//Choose one.

	if (cdrs_to_design.size() == 0){
		TR << "All CDRs fixed for low res graft designer or there are no CDRs in the set...." << std::endl;
		//Make sure our top designs has something in it for the antibody designer.
		check_for_top_designs(pose);
		return;
	}

	if (extend_native_cdrs_){
		extend_native_cdrs(pose, cdrs_to_design);
	}
	graft_mover_ = new CCDEndsGraftMover(ab_info_->get_CDR_start(cdrs_to_design[1], pose)-1, ab_info_->get_CDR_end(cdrs_to_design[1], pose)+1);

	//TR << "Superimposing ALL CDRs to graft onto target."<<std::endl;
	//for (core::Size i=1; i<=cdrs_to_design.size(); ++i){
		//CDRNameEnum cdr = cdrs_to_design[i];
		//core::Size start = ab_info_->get_CDR_start(cdr, pose)-1;
		//core::Size end = ab_info_->get_CDR_end(cdr, pose)+1;


		//for (core::Size pose_i = 1; pose_i <= cdr_set_[cdr].size(); ++pose_i){
		//	graft_mover_->set_insert_region(start, end);
		//	graft_mover_->set_piece(*(cdr_set_[cdr][pose_i]), overhang_, overhang_);
		//	graft_mover_->superimpose_overhangs_heavy(pose, false , false); //Superimposes first 4 atoms (Should be BB only then.)
		//}
	//}

	//Graft Settings.  Get optimal here!
	set_default_graft_settings();


	modeler_->set_scorefunction(scorefxn_);

	total_permutations_ = 1;
	for (core::Size i=1; i<=cdrs_to_design.size(); ++i){
		total_permutations_ *= cdr_set_[cdrs_to_design[i]].size();
	}
	TR<< "///// Total CDRs in set /////"<<std::endl;
	for(core::Size i = 1; i<=6; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		TR << "/// "<<ab_info_->get_CDR_name(cdr)<<" "<<cdr_set_[cdr].size()<<std::endl;
	}

	TR <<"Total possible CDR combinations: "<< total_permutations_ << std::endl;
	mc_ = new protocols::moves::MonteCarlo(pose, *scorefxn_, 1.0);

	core::Real native_score = (*scorefxn_)(pose);
	scorefxn_->show(pose);

	//if ((total_permutations_ <= graft_rounds_) && (use_deterministic_algorithm_ == true)  && (graft_rounds_ <= 10000)){
	//	run_deterministic_graft_algorithm(pose);
	//}
	//Temporary quick fix.
	if (cdrs_to_design.size() == 1 && use_deterministic_algorithm_){
		run_deterministic_graft_algorithm(pose, cdrs_to_design);
	}
	else{
		//run_random_graft_algorithm(pose, cdrs_to_design);
		run_basic_mc_algorithm(pose, cdrs_to_design);
	}

	TR << "Native Pose: " << native_score << std::endl;
	TR << "Final Pose: " << (*scorefxn_)(pose) << std::endl;
	for (core::Size i = 2; i<= top_scores_.size(); ++i){
		TR << "Top Ensemble " << i << " : " << (*scorefxn_)(*top_designs_[i]) << std::endl;
	}

	//Reinitialize ab_info clusters
	//ab_info_->setup_CDR_clusters(pose);

}

////////////////////////////////////////////// Boiler Plate ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string
AntibodyGraftDesignMover::get_name() const{
	return "AntibodyGraftDesignMover";
}


} //design
} //antibody
} //protocols
