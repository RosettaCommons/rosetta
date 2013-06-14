// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/design/AntibodyGraftDesigner.cc
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Project Includes
#include <protocols/antibody/design/AntibodyGraftDesigner.hh>
#include <protocols/antibody/CDRClusterEnum.hh>
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
#include <core/pose/PDBInfo.hh>

// Protocol Includes
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/grafting/AnchoredGraftMover.hh>
#include <protocols/grafting/util.hh>

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


static basic::Tracer TR("protocols.antibody.design.AntibodyGraftDesinger");
static numeric::random::RandomGenerator RG(1365421);

namespace protocols{
namespace antibody{
namespace design{
	using namespace protocols::antibody;
	using namespace protocols::grafting;
	using namespace core;
	using core::Size;
	
AntibodyGraftDesigner::AntibodyGraftDesigner(AntibodyInfoOP & ab_info):
	graft_mover_(NULL),
	scorefxn_(NULL)
{
	overhang_ = 3;
	ab_info_=ab_info;
	modeler_ = new AntibodyDesignModeler(ab_info_);
	if (ab_info_->get_Current_AntibodyNumberingScheme()!="Modified_AHO"){
		utility_exit_with_message("Antibody Design Protocol requires the Modified AHO numbering scheme");
	}
	set_defaults();
}

//AntibodyGraftDesigner::AntibodyGraftDesigner(AntibodyInfoOP & ab_info, std::string instruction_path) :
	//ab_info_(ab_info),
//	instruction_path_(instruction_path),
//	graft_mover_(NULL),
//	scorefxn_(NULL)
		
//{
//}

AntibodyGraftDesigner::~AntibodyGraftDesigner(){}

void
AntibodyGraftDesigner::set_defaults(){
	///Conservative defaults.  Defaults here are also read and set from a database file.  To allow run-time manipulations and testing.
	read_command_line_options();
	set_cdr_range(CDRNameEnum_start, CDRNameEnum_total, true);
	//set_cluster_range(CDRClusterEnum_start, CDRClusterEnum_stop, true);
	set_mintype_range(CDRNameEnum_start, CDRNameEnum_total, minimize);
	set_cdr_range_stay_in_native_cluster(CDRNameEnum_start, CDRNameEnum_total, false);
	set_cdr_range_stay_in_type(CDRNameEnum_start, CDRNameEnum_total, 1, true);
	set_cdr_cluster_centers_only_range(CDRNameEnum_start, CDRNameEnum_total, false);
	set_cdr_min_length_range(CDRNameEnum_start, CDRNameEnum_total, 1);
	set_cdr_max_length_range(CDRNameEnum_start, CDRNameEnum_total, 50); //Something absurd.
	set_keep_top_designs(10);
	read_instructions(instruction_path_);//Any settings in this file overwrite these defaults.  
}

void
AntibodyGraftDesigner::read_command_line_options(){
	set_graft_rounds(basic::options::option [basic::options::OptionKeys::antibody::design::graft_rounds]());
	instruction_path_ = basic::options::option [basic::options::OptionKeys::antibody::design::graft_instructions]();
}
void
AntibodyGraftDesigner::set_scorefunction(ScoreFunctionOP & scorefxn){
	scorefxn_=scorefxn->clone();
	if (scorefxn_->get_weight(core::scoring::dihedral_constraint)==0){
		scorefxn_->set_weight(core::scoring::dihedral_constraint, 1.0);
	}
}

void
AntibodyGraftDesigner::set_keep_top_designs(core::Size top_designs){
	num_top_designs_ = top_designs;
}

void
AntibodyGraftDesigner::set_graft_rounds(core::Size graft_rounds){
	graft_rounds_=graft_rounds;
}

void
AntibodyGraftDesigner::set_cluster(const CDRNameEnum cdr, const CDRClusterEnum cluster, const bool setting){
	if (setting==true){
		cdr_instructions_[cdr].leave_out_clusters.push_back(cluster);
	}
	else{
		cdr_instructions_[cdr].include_only_clusters.push_back(cluster);
	}
}

void
AntibodyGraftDesigner::set_cdr_set(CDRSet& cdr_set, CDRClusterMap & cdr_cluster_map, core::Size overhang){
	cdr_set_ = cdr_set;
	cdr_cluster_map_ = cdr_cluster_map;
	overhang_ = overhang;
}

void
AntibodyGraftDesigner::set_cluster_range(
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
AntibodyGraftDesigner::set_mintype(CDRNameEnum cdr_name, MinTypeEnum mintype){
	cdr_instructions_[cdr_name].mintype = mintype;
}

void
AntibodyGraftDesigner::set_mintype_range(CDRNameEnum cdr_start, CDRNameEnum cdr_end, MinTypeEnum mintype){
	for (core::SSize i=cdr_start; i<=cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		set_mintype(cdr, mintype);
	}
}

void
AntibodyGraftDesigner::set_cdr(const CDRNameEnum cdr, const bool setting){
	//cdrs_to_use_[cdr]=setting;
	cdr_instructions_[cdr].graft= setting;
}

void
AntibodyGraftDesigner::set_cdr_range(const CDRNameEnum cdr_start, CDRNameEnum cdr_end, const bool setting){
	for (core::SSize i=cdr_start; i<=cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		cdr_instructions_[cdr].graft = setting;
	}
}

void
AntibodyGraftDesigner::set_cdr_stay_in_native_cluster(const CDRNameEnum cdr, const bool setting){
	cdr_instructions_[cdr].stay_native_cluster = setting;
}

void
AntibodyGraftDesigner::set_cdr_range_stay_in_native_cluster(

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
AntibodyGraftDesigner::set_cdr_stay_in_type(const CDRNameEnum cdr, const core::Size type, const bool setting){
	cdr_instructions_[cdr].cluster_types[type] = setting;
}

void
AntibodyGraftDesigner::set_cdr_range_stay_in_type(

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
AntibodyGraftDesigner::set_cdr_min_length(CDRNameEnum const cdr, core::Size length){
	cdr_instructions_[cdr].min_length = length;
}

void
AntibodyGraftDesigner::set_cdr_min_length_range(const CDRNameEnum cdr_start, CDRNameEnum cdr_end, core::Size length){
	for (core::SSize i=cdr_start; i<=cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		set_cdr_min_length(cdr, length);
	}
}

void
AntibodyGraftDesigner::set_cdr_max_length(CDRNameEnum const cdr, core::Size length){
	cdr_instructions_[cdr].max_length = length;
}

void
AntibodyGraftDesigner::set_cdr_max_length_range(const CDRNameEnum cdr_start,const CDRNameEnum cdr_end, core::Size length){
	for (core::SSize i=cdr_start; i<=cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		set_cdr_max_length(cdr, length);
	}
}

void
AntibodyGraftDesigner::set_cdr_cluster_centers_only(const CDRNameEnum cdr, bool setting){
	cdr_instructions_[cdr].cluster_centers_only = setting; 
}

void
AntibodyGraftDesigner::set_cdr_cluster_centers_only_range(const CDRNameEnum cdr_start, CDRNameEnum cdr_end, bool setting){
	for (core::SSize i=cdr_start; i<=cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		set_cdr_cluster_centers_only(cdr, setting);
	}
}

void
AntibodyGraftDesigner::setup_native_clusters(pose::Pose & pose){
	if (! ab_info_->clusters_setup() ){
		ab_info_->setup_CDR_clusters(pose);
	}
	
	for (Size i=1; i<=CDRNameEnum_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		std::pair<CDRClusterEnum, core::Real> cluster_pair = ab_info_->get_CDR_cluster(cdr);
		native_clusters_[cdr] = cluster_pair.first;
	}
}

void
AntibodyGraftDesigner::initialize_cdr_set(){
	AntibodyDatabaseManagerOP manager = new AntibodyDatabaseManager();
	std::pair< CDRSet,  CDRClusterMap > result_pair = manager->load_cdrs_for_grafting(ab_info_, cdr_instructions_, pdbmap_);
	cdr_set_ = result_pair.first; cdr_cluster_map_ = result_pair.second;
	
}

void
AntibodyGraftDesigner::read_instructions(std::string instruction_path){
	DesignInstructionsParser parser = DesignInstructionsParser(ab_info_, instruction_path);
	parser.read_cdr_graft_instructions(cdr_instructions_);
}

void
AntibodyGraftDesigner::fix_pdb_info(pose::Pose pose, CDRNameEnum cdr, CDRClusterEnum cluster, core::Size original_start, core::Size original_pdb_end) {
	
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
		TR << "Set PDBInfo: "<< pose_num<<": "<< pdb_num << " "<<  chain <<std::endl;
		pose.pdb_info()->set_resinfo(pose_num, chain , pdb_num);
	}
	
	//Ending residues
	for (core::Size i = 1; i<=floor(result); ++i){
		core::Size pdb_num = pdb_num_end-i;
		core::Size pose_num = original_end-i;
		char chain = ab_info_->get_CDR_chain(cdr);
		TR << "Set PDBInfo: "<< pose_num<<": "<< pdb_num << " "<<  chain <<std::endl;
		pose.pdb_info()->set_resinfo(pose_num, chain , pdb_num);
	}
	
}

void
AntibodyGraftDesigner::graft_cdr(pose::Pose& pose, CDRNameEnum cdr, core::Size index){
	
	//Get CDR start/end in PDBInfo is now dynamic.  As long as pdb_info is not obsolete
	TR << "Grafting " << pdbmap_[cdr][index] << std::endl;
	TR <<  "Cluster " << ab_info_->get_cluster_name(cdr_cluster_map_[cdr][index]) << std::endl;
	core::Size start = ab_info_->get_CDR_start(cdr, pose)-1;
	core::Size end = ab_info_->get_CDR_end(cdr, pose)+1;
	//core::Size original_pdb_end = pose.pdb_info()->number(end);
	
	graft_mover_->set_insert_region(start, end);
	//(cdr_set_[cdr][index])->dump_pdb("piece.pdb");
	pose::Pose piece = pose::Pose(*(cdr_set_[cdr][index]));
	graft_mover_->set_piece(piece, overhang_, overhang_);
	if (cdr==h3){
		//graft_mover_->set_scaffold_flexibility(3, 3);	
		graft_mover_->set_cycles(40);
	}
	graft_mover_->superimpose_overhangs_heavy(pose, false , true); //Superimposes first 4 atoms (Should be BB only then.)
	graft_mover_->apply(pose);
	
	
	//fix_pdb_info(pose, cdr, cdr_cluster_map_[cdr][index], start, original_pdb_end);//Needs to run after every graft before constraints can be added properly
	pose.pdb_info()->copy(*((cdr_set_[cdr][index])->pdb_info()), 1 + overhang_, (cdr_set_[cdr][index])->total_residue()-overhang_, start+1);
	pose.pdb_info()->obsolete( false );
	modeler_->set_cdr_only(cdr, true);
	
	
	///Minimize sidechains or cdr
	if (cdr_instructions_[cdr].mintype == relax){
		graft_mover_->repack_connection_and_residues_in_movemap(pose, scorefxn_);
		protocols::antibody::set_harmonic_constraint(ab_info_, pose, cdr_cluster_map_[cdr][index]);
		modeler_->relax_cdrs(pose, scorefxn_, false);//Change this if need be.
		//Include close antigen sidechains here?
	}
	else if (cdr_instructions_[cdr].mintype == centroid_relax){
		graft_mover_->repack_connection_and_residues_in_movemap(pose, scorefxn_);
		protocols::antibody::set_harmonic_constraint(ab_info_, pose, cdr_cluster_map_[cdr][index]);
		modeler_->relax_cdrs(pose, scorefxn_, true);//Change this if need be.
	}
	else if (cdr_instructions_[cdr].mintype == repack){
		graft_mover_->repack_connection_and_residues_in_movemap_and_piece(pose, scorefxn_);
	}
	else if (cdr_instructions_[cdr].mintype == minimize){
		graft_mover_->repack_connection_and_residues_in_movemap_and_piece(pose, scorefxn_);
		protocols::antibody::set_harmonic_constraint(ab_info_, pose, cdr_cluster_map_[cdr][index]);
		modeler_->minimize_cdrs(pose, scorefxn_);
		
	}
	else {
		graft_mover_->repack_connection_and_residues_in_movemap(pose, scorefxn_);
	}
	
	check_for_top_designs(pose);
	
}

void
AntibodyGraftDesigner::check_for_top_designs(pose::Pose & pose){
	
	//Can be refactored to use utility::TopScoreSelector
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

void
AntibodyGraftDesigner::run_random_graft_algorithm(pose::Pose & pose, vector1<CDRNameEnum> & cdrs_to_design){
	TR <<"Running default graft algorithm"<<std::endl;
	
	
	//First version of algorithm.  Will eventually use DEE/Tree decomposition to sample. 
	
	//Map and vector place it's elements by default in the Heap.  So, this should be OK.
	std::map< vector1< Size >, Size  >used_indexes;
	vector1< vector1< Size > >all_indexes;
	 
	//Size value represents the number of grafts it took to accomplish this combination.
	top_scores_.push_back((*scorefxn_)(pose));
	top_designs_.push_back(new Pose());
	*(top_designs_[1]) = pose;
	
	//Cannot create index list and choose from that and pop it.  Too many possibilities.  Not enough RAM. (6x10^15 numbers if maximal sampling)
	//We also do not want to graft more then one CDR at a time if we can avoid it.  
	bool keep_picking;;
	Size total_trials = 300;
	
	vector1< Size >current_index (6, 0); //These represent the CDRs we have on the current pose
	vector1<Size>regraft_index;
	
	if (graft_rounds_ >= 1000000){
		utility_exit_with_message("Grafting was not designed for this many rounds.  Please use less rounds");
	}
	//Get indexes to choose from since memory will be low.
	if (total_permutations_ <= 10000){
		TR.Debug<<"Begin loading Permutation indexes"<<std::endl;
		vector1<core::Size> total_cdr_set(6, 0);
		for (Size i=1; i<=cdrs_to_design.size(); ++i){
			CDRNameEnum cdr = cdrs_to_design[i];
			total_cdr_set[i]=cdr_set_[cdr].size();
			get_all_graft_permutations(total_cdr_set, all_indexes, current_index, 1);
		}
		TR.Debug<<"All Permutation indexes loaded"<<std::endl;
		
	}
	for (Size i = 1; i<=graft_rounds_; ++i){
		TR << "Graft round: " << i <<std::endl;
		Size cdr_index;
		Size pose_index;
		Size n_trials = 0;
		keep_picking=true;
		
		while(keep_picking){
			
			//Pick random CDR and random CDR to graft.
			cdr_index = RG.random_range(1, cdrs_to_design.size());
			CDRNameEnum cdr = cdrs_to_design[cdr_index];
			pose_index = RG.random_range(1, cdr_set_[cdr].size());
			
			//Create new index and test if it has been seen
			vector1< Size> new_index(current_index);
			new_index[cdr_index] = pose_index;
			std::map< vector1< Size >, Size>::const_iterator iter(used_indexes.find(new_index));
			
			if (iter == used_indexes.end()){
				//->Combination found ok.
				keep_picking=false;
			} 
			else if(n_trials == total_trials){
				//->Would take many trials to find one that hasn't been done.  So do something else.  Find one that hasn't been used.
				//Try to only graft ONE CDR if possible.  We don't want to graft many if we can help it.
				TR << "Perumutations Left: " << total_permutations_-used_indexes.size() << "Rounds Left: "<< graft_rounds_-i << std::endl;
				for (Size z = 1; z<=cdrs_to_design.size(); ++z){
					CDRNameEnum cdr = cdrs_to_design[z];
					vector1< Size > temp_index (current_index);
					for (Size j=1; j<=cdr_set_[cdr].size(); ++j){
						temp_index[i]=j;
						std::map< vector1< Size >, Size>::const_iterator iter(used_indexes.find(temp_index));
						//Check if we have used the combination.
						if (iter==used_indexes.end()){
							keep_picking=false;
							break;
						} 
					}
					if (! keep_picking){break;}
				}
				
				
				while(keep_picking){
					if(all_indexes.size() != 0){
						core::Size ind = RG.random_range(1, all_indexes.size());
						regraft_index = all_indexes[ind];
						keep_picking=false;
					}
					else {
						//Choose random index of each CDR instead of just 1.  Test if it has been used. 
						regraft_index =current_index;
						cdr_index = RG.random_range(1, cdrs_to_design.size());
						cdr = cdrs_to_design[cdr_index];
						pose_index = RG.random_range(1, cdr_set_[cdr].size());
						regraft_index[cdr_index]=pose_index;
						std::map< vector1< Size >, Size>::const_iterator iter(used_indexes.find(regraft_index));
						if (iter==used_indexes.end()){
							keep_picking=false;
						}
					}
				}
			} //ntrials = total_trials
			else{
				n_trials += 1;//Keep trying
			}
		}
		if (regraft_index.size() != 0){
			//Here, we need to graft more than one CDR
			Size grafted_cdrs = 0;
			for (core::Size x = 1; x<=cdrs_to_design.size(); ++i){
				if(current_index[x]!= regraft_index[x]){
					CDRNameEnum cdr = cdrs_to_design[x];
					graft_cdr(pose, cdr, regraft_index[x]);
					grafted_cdrs+=1;
				}
			}
			//Apply boltzmann criterion.  Set current index if accepted.
			if (mc_->boltzmann(pose)){
				used_indexes[regraft_index]=grafted_cdrs;
				current_index = regraft_index;
			}
			
			regraft_index.clear();
		}
		else{
			
			//Normal way - Only need to graft 1 CDR
			CDRNameEnum cdr = cdrs_to_design[cdr_index];
			graft_cdr(pose, cdr, pose_index);
			//Apply boltzmann criterion.  Set current index if accepted.
			if (mc_->boltzmann(pose)){
				vector1< Size> new_index(current_index);
				new_index[cdr_index] = pose_index;
				used_indexes[new_index]=1;
				current_index =new_index;
			}
		}
		//Take care of all_indexes if present.
		if (all_indexes.size() != 0 ){
			vector1< vector1< core::Size > >::Iterator iter = std::find(all_indexes.begin(), all_indexes.end(), current_index);
			if (iter !=  all_indexes.end()){
				all_indexes.erase(iter);
			}
		}
	}
	mc_->recover_low(pose);
}


void
AntibodyGraftDesigner::run_deterministic_graft_algorithm(pose::Pose & pose, vector1<CDRNameEnum> & cdrs_to_design, core::Size recurse_num){
	//Note:  Not feasible with >= 4 CDRs to try at each position.
	TR<< "Running deterministic graft algorithm -  Trying every possible CDR combination"<<std::endl;
	
	if (recurse_num > cdrs_to_design.size()){
		return;
	} 
	else if (recurse_num==0){

		top_scores_.push_back((*scorefxn_)(pose));
		top_designs_.push_back(new Pose());
		*(top_designs_[1]) = pose;
		run_deterministic_graft_algorithm(pose, cdrs_to_design, 1);
		mc_->recover_low(pose);	
		
	} else{
		//Recursive due to variable number of nested For Loops
		for (core::Size i=1; i<=cdr_set_[cdrs_to_design[recurse_num]].size(); ++i){
			CDRNameEnum cdr = cdrs_to_design[recurse_num];
			
			graft_cdr(pose, cdr, i);
			mc_->eval_lowest_score_pose(pose, false, true);
			if (recurse_num == cdrs_to_design.size()){
				//We are all the way in
				continue;
			}
			run_deterministic_graft_algorithm(pose, cdrs_to_design, recurse_num+1);
		}
	}
}

void
AntibodyGraftDesigner::apply(pose::Pose & pose){
	
	//Create Instances.  Make sure everything is ready to begin.
	if (!scorefxn_){
		scorefxn_ = core::scoring::getScoreFunction(true);
		if (scorefxn_->get_weight(core::scoring::dihedral_constraint)==0.0){
			scorefxn_->set_weight(core::scoring::dihedral_constraint, 1.0);
		}
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
		TR << "//////// " << ab_info_->get_CDR_Name(cdr) << "///////////////////////////////////////////////"<< std::endl;
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
	
	//List of CDRs to graft.
	vector1< CDRNameEnum > cdrs_to_design;
	
	for (core::Size i=1; i<=CDRNameEnum_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		if (cdr_instructions_[cdr].graft && cdr_set_[cdr].size() != 0){
			cdrs_to_design.push_back(cdr);
		}
	}
	//Choose one.
	
	//Setup Graft + Superimpose ALL CDRs first.  Instead of doing it everytime we need to graft.
	
	graft_mover_ = new AnchoredGraftMover(ab_info_->get_CDR_start(cdrs_to_design[1], pose)-1, ab_info_->get_CDR_end(cdrs_to_design[1], pose)+1);
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
	graft_mover_->set_cycles(15);
	graft_mover_->set_use_smooth_centroid_settings(true);
	graft_mover_->set_use_double_loop_double_CCD_arms(true);
	
	total_permutations_ = 1;
	for (core::Size i=1; i<=cdrs_to_design.size(); ++i){
		total_permutations_ *= cdr_set_[cdrs_to_design[i]].size();
	}
	TR<< "///// Total CDRs /////"<<std::endl;
	for(core::Size i = 1; i<=6; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		TR << "/// "<<ab_info_->get_CDR_Name(cdr)<<" "<<cdr_set_[cdr].size()<<std::endl;
	}
	
	TR <<"Total possible CDR combinations: "<< total_permutations_ << std::endl;
	mc_ = new protocols::moves::MonteCarlo(pose, *scorefxn_, 1.0);
	
	core::Real native_score = (*scorefxn_)(pose);
	if (total_permutations_<=graft_rounds_){
		run_deterministic_graft_algorithm(pose, cdrs_to_design, 0);
	}
	else{
		run_random_graft_algorithm(pose, cdrs_to_design);
	}
	TR << "Native Pose: " << native_score << std::endl;
	TR << "Final Pose: " << (*scorefxn_)(pose) << std::endl;
	for (core::Size i = 2; i<= top_scores_.size(); ++i){
		TR << "Top Ensemble " << i << " : " << (*scorefxn_)(*top_designs_[i]) << std::endl;
	}
	
	//Reinitialize ab_info clusters
	ab_info_->setup_CDR_clusters(pose);
	
}

////////////////////////////////////////////// Boiler Plate ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string
AntibodyGraftDesigner::get_name() const{
	return "AntibodyGraftDesigner";
}


} //design
} //antibody
} //protocols
