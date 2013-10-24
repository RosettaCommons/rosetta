// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody_design/AntibodyDesignModeler.cc
/// @brief Handles modeling of the antibody.  Before and after design
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

//Antibody Headers
#include <protocols/antibody/design/AntibodyDesignModeler.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/CDRClusterEnum.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>

// Core Headers
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreType.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>

//Protocol Headers
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/CentroidRelax.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/docking/DockingLowRes.hh>
#include <protocols/docking/DockMCMProtocol.hh>
#include <protocols/docking/util.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <protocols/toolbox/task_operations/PreventChainFromRepackingOperation.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
//#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/file/FileName.hh>
#include <sstream>

static basic::Tracer TR("antibody.design.AntibodyDesignModeler");
namespace protocols {
namespace antibody {
namespace design {
            
            
using namespace core::kinematics;         
using namespace core::pose;    
using namespace core::scoring;           
using namespace protocols::simple_moves;
using namespace protocols::toolbox::task_operations;
using core::Size;      
using std::string;
            
        
AntibodyDesignModeler::AntibodyDesignModeler(AntibodyInfoOP ab_info){
	ab_info_ = ab_info;
	cdrs_.resize(CDRNameEnum_total);
	set_cdr_range(CDRNameEnum_start, CDRNameEnum_total, true);
	read_command_line_options();
	
	scorefxn_ = getScoreFunction();
	docking_scorefxn_high_ = ScoreFunctionFactory::create_score_function( "docking", "docking_min" ) ;
	cmd_line_operation_ = new core::pack::task::operation::InitializeFromCommandline();
	restrict_design_operation_ = new core::pack::task::operation::RestrictToRepacking();
	loops_operation_ = new RestrictToLoopsAndNeighbors();
	
	tf_ = new core::pack::task::TaskFactory();
	antigen_interface_tf_ = new core::pack::task::TaskFactory();
	antigen_interface_tf_->push_back(cmd_line_operation_);
	antigen_interface_tf_->push_back( restrict_design_operation_);
	antigen_interface_tf_->push_back(new RestrictToInterface(1, interface_dis_/*Distance*/));
	
	
}
        
        
AntibodyDesignModeler::~AntibodyDesignModeler(){}
        
void
AntibodyDesignModeler::read_command_line_options(){
	set_interface_detection_dis(basic::options::option [basic::options::OptionKeys::antibody::design::interface_dis]());
	set_neighbor_detection_dis(basic::options::option [basic::options::OptionKeys::antibody::design::neighbor_dis]());
	set_ab_dock_chains(basic::options::option [basic::options::OptionKeys::antibody::design::ab_dock_chains]());
}

protocols::loops::LoopsOP
AntibodyDesignModeler::get_cdr_loops(Pose& pose) const {
	
	protocols::loops::LoopsOP cdr_loops = new protocols::loops::Loops;
	
	for (Size i = 1; i <= CDRNameEnum_total; ++i){
		if (cdrs_[i]){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			Size start = ab_info_->get_CDR_start(cdr, pose);
			Size stop =  ab_info_->get_CDR_end(cdr, pose);
			Size cutpoint = (stop-start)/2+start;
			protocols::loops::Loop cdr_loop = protocols::loops::Loop(start, stop, cutpoint);
			cdr_loops->add_loop(cdr_loop);
		}
	}
	return cdr_loops;
}

void
AntibodyDesignModeler::apply_LH_A_foldtree(core::pose::Pose & pose) const {
	vector1< char > antigen_chains = ab_info_->get_antigen_chains();
	
	std::string antigen(antigen_chains.begin(), antigen_chains.end());
	std::string dock_chains = "LH_"+antigen;
	
	vector1< int > movable_jumps(1, 1);
	protocols::docking::setup_foldtree(pose, dock_chains, movable_jumps);
}

void
AntibodyDesignModeler::set_cdr_range(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, bool setting) {
	for (core::SSize i=cdr_start; i<=cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		set_cdr(cdr, setting);
	}
}

void
AntibodyDesignModeler::set_cdr_only(CDRNameEnum cdr, bool setting){
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
AntibodyDesignModeler::set_cdr(CDRNameEnum const cdr, bool setting){
	cdrs_[cdr]=setting;
}

void
AntibodyDesignModeler::set_ab_dock_chains(std::string ab_dock_chains){
	ab_dock_chains_ = ab_dock_chains;
}

std::string
AntibodyDesignModeler::get_ab_dock_chains(){
	return ab_dock_chains_;
}

void
AntibodyDesignModeler::set_scorefunction(ScoreFunctionCOP scorefxn){
	scorefxn_ = scorefxn->clone();
}

void
AntibodyDesignModeler::set_interface_detection_dis(core::Real interface_distance){
	interface_dis_ = interface_distance;
}

void
AntibodyDesignModeler::set_neighbor_detection_dis(core::Real neighbor_distance){
	neighbor_dis_ = neighbor_distance;
}

void
AntibodyDesignModeler::relax_cdrs(Pose & pose, bool centroid_mode, bool starting_coordinate_constraints /*false*/, bool min_interface /* false */, std::string dock_chains /* L_H*/ ) const {

	TR << "Relaxing CDR(s)" << std::endl;
	ScoreFunctionOP local_scorefxn = scorefxn_->clone();
	
	MoveMapOP mm = new MoveMap();
	local_scorefxn->set_weight(chainbreak, 100.00); //First PyRosetta Toolkit lesson.
	TR << "start: "<<(*local_scorefxn)(pose) << std::endl;
	core::kinematics::FoldTree original_ft = pose.fold_tree();
	if (min_interface){
		TR << "Including RB interface min" <<std::endl;
		vector1< int > movable_jumps(1, 1);
		protocols::docking::setup_foldtree(pose, dock_chains, movable_jumps);
	}
	for (core::Size i=1; i<=CDRNameEnum_total; ++i){
		if (cdrs_[i]){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			ab_info_->add_CDR_to_MoveMap(pose, mm, cdr);
		}
	}
	protocols::loops::LoopsOP cdr_loops = get_cdr_loops(pose);
	loops_operation_->set_design_loop(false);
	loops_operation_->set_include_neighbors(false);
	loops_operation_->set_loops(cdr_loops);
	loops_operation_->set_cutoff_distance(neighbor_dis_);
	
	tf_->clear();
	tf_->push_back(cmd_line_operation_ );
	tf_->push_back( restrict_design_operation_ );
	tf_->push_back(loops_operation_);
	
	if (centroid_mode){
		protocols::relax::CentroidRelaxOP rel = new protocols::relax::CentroidRelax(mm);
		if (starting_coordinate_constraints){
			rel->constrain_relax_to_start_coords(true);
		}
		rel->apply(pose);
	}
	else{
		
		protocols::relax::FastRelaxOP rel = new protocols::relax::FastRelax(local_scorefxn);//Stack construction failed to compile.
		rel->set_movemap(mm);
		rel->set_task_factory(tf_);
		if (starting_coordinate_constraints){
			rel->constrain_relax_to_start_coords(true);
		}
		rel->apply(pose);
		
	}
	TR << "end: "<<(*local_scorefxn)(pose) << std::endl;
	if (starting_coordinate_constraints){
		core::scoring::constraints::remove_constraints_of_type(pose, "CoordinateConstraint");
	}
	pose.fold_tree(original_ft);
}

void
AntibodyDesignModeler::relax_cdrs_and_neighbor_sc(Pose & pose, bool starting_coordinate_constraints /*false*/, bool min_interface /* false */, std::string dock_chains /* L_H*/ ) const {
	
	TR << "Relaxing CDR(s) and Neighbor Sidechains" << std::endl;
	ScoreFunctionOP local_scorefxn = scorefxn_->clone();
	local_scorefxn->set_weight(chainbreak, 100.00); //First PyRosetta Toolkit lesson.
	
	core::kinematics::FoldTree original_ft = pose.fold_tree();
	if (min_interface){
		TR << "Including RB interface min" <<std::endl;
		vector1< int > movable_jumps(1, 1);
		protocols::docking::setup_foldtree(pose, dock_chains, movable_jumps);
	}
	
	TR << "start: "<<(*local_scorefxn)(pose) << std::endl;
	MoveMapOP mm = new MoveMap();
	for (core::Size i=1; i<=CDRNameEnum_total; ++i){
		if (cdrs_[i]){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			ab_info_->add_CDR_to_MoveMap(pose, mm, cdr, false, true, neighbor_dis_);
		}
	}
	protocols::loops::LoopsOP cdr_loops = get_cdr_loops(pose);
	loops_operation_->set_design_loop(false);
	loops_operation_->set_include_neighbors(true);
	loops_operation_->set_loops(cdr_loops);
	loops_operation_->set_cutoff_distance(neighbor_dis_);
	
	tf_->clear();
	tf_->push_back(cmd_line_operation_ );
	tf_->push_back( restrict_design_operation_ );
	tf_->push_back(loops_operation_);
	
	protocols::relax::FastRelaxOP rel = new protocols::relax::FastRelax(local_scorefxn);//Stack construction failed to compile.
	rel->set_movemap(mm);
	rel->set_task_factory(tf_);
	if (starting_coordinate_constraints){
		rel->constrain_relax_to_start_coords(true);
	}
	rel->apply(pose);
	TR << "end: "<<(*local_scorefxn)(pose) << std::endl;
	
	if (starting_coordinate_constraints){
		core::scoring::constraints::remove_constraints_of_type(pose, "CoordinateConstraint");
	}
	
	pose.fold_tree(original_ft);
}

void
AntibodyDesignModeler::relax_interface(Pose & pose, std::string dock_chains, bool min_interface_sc /* true */) const {
	
	TR << "Relaxing Interface: " << dock_chains << std::endl;
	ScoreFunctionOP local_scorefxn = scorefxn_->clone();
	vector1< int > movable_jumps(1, 1);
	core::kinematics::FoldTree original_ft = pose.fold_tree();
	protocols::docking::setup_foldtree(pose, dock_chains, movable_jumps);
	
	MoveMapOP mm = new MoveMap();
	mm->set_jump(true); //Should be only one jump
	core::pack::task::TaskFactoryOP tf = antigen_interface_tf_->clone();
	core::pack::task::PackerTaskOP task = antigen_interface_tf_->create_task_and_apply_taskoperations(pose);
	protocols::relax::FastRelaxOP rel = new protocols::relax::FastRelax(local_scorefxn);
	
	//Add residues to movemap, so that they minimize during minimization, while repacking will happen via tf. Maybe should not have any chi packing.
	if (min_interface_sc){
		for (core::Size i = 1; i<=pose.total_residue(); i++){
			if (task->pack_residue(i)){
				mm->set_chi(i, true);
			}
		}
	}
	
	rel->set_movemap(mm);
	rel->set_task_factory(tf);
	rel->apply(pose);
	pose.fold_tree(original_ft);
}

void
AntibodyDesignModeler::minimize_interface(Pose& pose, std::string dock_chains, bool min_interface_sc /* true */) const {
	
	TR << "Minimizing Interface: " << dock_chains << std::endl;
	vector1< int > movable_jumps(1, 1);
	core::kinematics::FoldTree original_ft = pose.fold_tree();
	protocols::docking::setup_foldtree(pose, dock_chains, movable_jumps);
	
	MoveMapOP mm = new MoveMap();
	mm->set_jump(true); //Should be only one jump
	
	core::pack::task::PackerTaskOP task = antigen_interface_tf_->create_task_and_apply_taskoperations(pose);
	if (min_interface_sc){
		for (core::Size i = 1; i<=pose.total_residue(); i++){
			if (task->pack_residue(i)){
				mm->set_chi(i, true);
			}
		}
	}
	ScoreFunctionOP local_scorefxn = scorefxn_->clone();
	protocols::simple_moves::MinMover min_mover(mm, local_scorefxn, "dfpmin_armijo_nonmonotone", 0.01, false /*use_nblist*/ );
	min_mover.apply(pose);
	
	
	pose.fold_tree(original_ft);
}

void
AntibodyDesignModeler::minimize_cdrs(Pose & pose, bool min_interface /* false */, std::string dock_chains /* L_H*/ ) const {
	
	TR << "Minimizing CDR(s)" << std::endl;
	MoveMapOP mm = new MoveMap();
	ScoreFunctionOP local_scorefxn = scorefxn_->clone();
	local_scorefxn->set_weight(chainbreak, 100.00); 
	TR << "start: "<<(*local_scorefxn)(pose) << std::endl;
	
	core::kinematics::FoldTree original_ft = pose.fold_tree();
	if (min_interface){
		TR << "Including RB interface min" <<std::endl;
		vector1< int > movable_jumps(1, 1);
		protocols::docking::setup_foldtree(pose, dock_chains, movable_jumps);
	}
	for (core::Size i = 1; i<=CDRNameEnum_total; ++i){
		if (cdrs_[i]){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			ab_info_->add_CDR_to_MoveMap(pose, mm, cdr);
		}
	}

	protocols::simple_moves::MinMover min_mover(mm, local_scorefxn, "dfpmin_armijo_nonmonotone", 0.01, false /*use_nblist*/ );
	min_mover.apply(pose);
	pose.fold_tree(original_ft);
	TR << "end: "<<(*local_scorefxn)(pose) << std::endl;
}

void
AntibodyDesignModeler::minimize_cdrs_and_neighbor_sc(Pose& pose, bool min_interface /* false */, std::string dock_chains /* L_H*/ ) const {
	
	TR << "Minimizing CDR(s) and neighbor Sidechains" << std::endl;
	MoveMapOP mm = new MoveMap();
	ScoreFunctionOP local_scorefxn = scorefxn_->clone();
	local_scorefxn->set_weight(chainbreak, 100.00);
	TR << "start: "<<(*local_scorefxn)(pose) << std::endl;
	
	core::kinematics::FoldTree original_ft = pose.fold_tree();
	if (min_interface){
		TR << "Including RB interface min" <<std::endl;
		vector1< int > movable_jumps(1, 1);
		protocols::docking::setup_foldtree(pose, dock_chains, movable_jumps);
	}
	for (core::Size i = 1; i<=CDRNameEnum_total; ++i){
		if (cdrs_[i]){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			ab_info_->add_CDR_to_MoveMap(pose, mm, cdr, false, true, neighbor_dis_);
		}
	}

	protocols::simple_moves::MinMover min_mover(mm, local_scorefxn, "dfpmin_armijo_nonmonotone", 0.01, false /*use_nblist*/ );
	min_mover.apply(pose);
	pose.fold_tree(original_ft);
	TR << "end: "<<(*local_scorefxn)(pose) << std::endl;
}

void
AntibodyDesignModeler::repack_antigen_ab_interface(Pose& pose) const {
	vector1< char > antigen_chains = ab_info_->get_antigen_chains();
	if (antigen_chains.size() == 0){
		TR <<" Antigen not present to repack interface" << std::endl;
		return;
	}
	
	TR << "Repacking antigen-antibody interface" <<std::endl;
	core::kinematics::FoldTree original_ft = pose.fold_tree();
	apply_LH_A_foldtree(pose);
	TR << "start: "<<(*scorefxn_)(pose) << std::endl;
	core::pack::task::PackerTaskOP task = antigen_interface_tf_->create_task_and_apply_taskoperations(pose);
	
	PackRotamersMover packer = PackRotamersMover(scorefxn_, task);
	packer.apply(pose);
	TR << "end: "<<(*scorefxn_)(pose) << std::endl;
	pose.fold_tree(original_ft);
}

void
AntibodyDesignModeler::repack_antigen_interface(Pose & pose) const {
	vector1< char > antigen_chains = ab_info_->get_antigen_chains();
	if (antigen_chains.size() == 0){
		TR <<" Antigen not present to repack interface" << std::endl;
		return;
	}
	
	TR << "Repacking antigen part of interface" <<std::endl;
	core::kinematics::FoldTree original_ft = pose.fold_tree();
	apply_LH_A_foldtree(pose);
	TR << "start: "<<(*scorefxn_)(pose) << std::endl;
	core::pack::task::TaskFactoryOP tf = antigen_interface_tf_->clone();
	Size L_chain  = core::pose::get_chain_id_from_chain('L', pose);
	Size H_chain = core::pose::get_chain_id_from_chain('H', pose);
	tf->push_back(new PreventChainFromRepackingOperation(L_chain));
	tf->push_back(new PreventChainFromRepackingOperation(H_chain));
	
	core::pack::task::PackerTaskOP task = tf->create_task_and_apply_taskoperations(pose);
	
	PackRotamersMover packer = PackRotamersMover(scorefxn_, task);
	packer.apply(pose);
	TR << "end: "<<(*scorefxn_)(pose) << std::endl;
	pose.fold_tree(original_ft);
}

void
AntibodyDesignModeler::repack_antibody_interface(Pose & pose) const {
	vector1< char > antigen_chains = ab_info_->get_antigen_chains();
	if (antigen_chains.size() == 0){
		TR <<" Antigen not present to detect antibody interface" << std::endl;
		return;
	}
	
	TR << "Repacking antibody part of interface" << std::endl;
	core::kinematics::FoldTree original_ft = pose.fold_tree();
	apply_LH_A_foldtree(pose);
	TR << "start: "<<(*scorefxn_)(pose) << std::endl;
	core::pack::task::TaskFactoryOP tf = antigen_interface_tf_->clone();
	
	for (Size i = 1; i <= antigen_chains.size(); ++i){
		Size chain = core::pose::get_chain_id_from_chain(antigen_chains[i], pose);
		tf->push_back(new PreventChainFromRepackingOperation(chain));
	}
	core::pack::task::PackerTaskOP task = tf->create_task_and_apply_taskoperations(pose);
	
	PackRotamersMover packer = PackRotamersMover(scorefxn_, task);
	packer.apply(pose);
	TR << "end: "<<(*scorefxn_)(pose) << std::endl;
	pose.fold_tree(original_ft);
}

void
AntibodyDesignModeler::repack_CDRs(Pose& pose) {
	using namespace protocols::toolbox::task_operations;
	
	TR << "Repacking CDR(s)" << std::endl;
	protocols::loops::LoopsOP cdr_loops = get_cdr_loops(pose);
	
	TR << "start: "<<(*scorefxn_)(pose) << std::endl;
	loops_operation_->set_design_loop(false);
	loops_operation_->set_include_neighbors(false);
	loops_operation_->set_loops(cdr_loops);

	tf_->clear();
	tf_->push_back(cmd_line_operation_ );
	tf_->push_back( restrict_design_operation_ );
	tf_->push_back(loops_operation_);
	
	core::pack::task::PackerTaskOP task = tf_->create_task_and_apply_taskoperations(pose);
	
	PackRotamersMover packer = PackRotamersMover(scorefxn_, task);
	packer.apply(pose);
	TR << "end: "<<(*scorefxn_)(pose) << std::endl;
	
}
 
void
AntibodyDesignModeler::repack_CDRs_and_neighbors(Pose& pose){
	using namespace protocols::toolbox::task_operations;
	
	TR << "Repacking CDR(s) and Neighbors" <<std::endl;
	TR << "start: "<<(*scorefxn_)(pose) << std::endl;
	protocols::loops::LoopsOP cdr_loops = get_cdr_loops(pose);
	loops_operation_->set_design_loop(false);
	loops_operation_->set_include_neighbors(true);
	loops_operation_->set_loops(cdr_loops);
	loops_operation_->set_cutoff_distance(neighbor_dis_);
	
	tf_->clear();
	tf_->push_back(cmd_line_operation_ );
	tf_->push_back( restrict_design_operation_ );
	tf_->push_back(loops_operation_);
	
	core::pack::task::PackerTaskOP task = tf_->create_task_and_apply_taskoperations(pose);
	
	PackRotamersMover packer = PackRotamersMover(scorefxn_, task);
	packer.apply(pose);
	TR << "end: "<<(*scorefxn_)(pose) << std::endl;
}

void
AntibodyDesignModeler::dock_LH_A_low_res(Pose& pose, bool pack_interface /*false*/) const {

	vector1< char > antigen_chains = ab_info_->get_antigen_chains();
	if (antigen_chains.size() == 0){
		TR <<" Antigen not present to dock" << std::endl;
		return;
	}
	
	
	std::string antigen(antigen_chains.begin(), antigen_chains.end());
	std::string dock_chains = ab_dock_chains_+"_"+antigen;
	TR << "Docking " <<dock_chains << " low res " << std::endl;
	TR.Debug << "dock_chains: "<< dock_chains << std::endl;
	
	dock_low_res(pose, dock_chains, pack_interface);
}

void
AntibodyDesignModeler::dock_LH_A_high_res(Pose & pose, int first_cycle /* 4 */, int second_cycle /* 45 */) const {

	vector1< char > antigen_chains = ab_info_->get_antigen_chains();
	if (antigen_chains.size() == 0){
		TR <<" Antigen not present to dock" << std::endl;
		return;
	}
	
	std::string antigen(antigen_chains.begin(), antigen_chains.end());
	std::string dock_chains = ab_dock_chains_+"_"+antigen;
	TR << "Docking " <<dock_chains << " high res " << std::endl;
	dock_high_res(pose, dock_chains, first_cycle, second_cycle);
		
}

void
AntibodyDesignModeler::dock_A_LH_low_res(Pose& pose, bool pack_interface /*false*/) const {
	vector1< char > antigen_chains = ab_info_->get_antigen_chains();
	if (antigen_chains.size() == 0){
		TR <<" Antigen not present to dock" << std::endl;
		return;
	}
	
	std::string antigen(antigen_chains.begin(), antigen_chains.end());
	std::string dock_chains = antigen+"_"+ab_dock_chains_;
	TR << "Docking " <<dock_chains << " low res " << std::endl;
	dock_low_res(pose, dock_chains, pack_interface);
}

void
AntibodyDesignModeler::dock_A_LH_high_res(Pose & pose, int first_cycle /* 4 */, int second_cycle /* 45 */) const {

	vector1< char > antigen_chains = ab_info_->get_antigen_chains();
	if (antigen_chains.size() == 0){
		TR <<" Antigen not present to dock" << std::endl;
		return;
	}
	
	std::string antigen(antigen_chains.begin(), antigen_chains.end());
	std::string dock_chains = antigen+"_"+ab_dock_chains_;
	TR << "Docking " <<dock_chains << " high res " << std::endl;
	dock_high_res(pose, dock_chains, first_cycle, second_cycle);
		
}

void
AntibodyDesignModeler::dock_L_H_low_res(Pose& pose, bool pack_interface /*false*/) const {

	TR << "Docking L-H low res" << std::endl;
	std::string dock_chains = "L_H";
	
	dock_low_res(pose, dock_chains, pack_interface);
}

void
AntibodyDesignModeler::dock_L_H_high_res(Pose & pose, int first_cycle /* 4 */, int second_cycle /* 45 */) const {

	TR << "Docking L-H high res " << std::endl;
	std::string dock_chains = "L_H";
	
	dock_high_res(pose, dock_chains, first_cycle, second_cycle);
		
}

void
AntibodyDesignModeler::dock_H_L_low_res(Pose& pose, bool pack_interface /*false*/) const {
	
	TR << "Docking H-L low res" << std::endl;
	std::string dock_chains = "H_L";
	
	dock_low_res(pose, dock_chains, pack_interface);
}

void
AntibodyDesignModeler::dock_H_L_high_res(Pose & pose, int first_cycle /* 4 */, int second_cycle /* 45 */) const {

	TR << "Docking H-L high res" << std::endl;
	std::string dock_chains = "H_L";
	
	dock_high_res(pose, dock_chains, first_cycle, second_cycle);
		
}

void
AntibodyDesignModeler::dock_high_res(Pose & pose, std::string dock_chains, int first_cycle /* 4 */, int second_cycle /* 45 */) const {

	
	vector1< int > movable_jumps(1, 1);
	core::kinematics::FoldTree original_ft = pose.fold_tree();
	protocols::docking::setup_foldtree(pose, dock_chains, movable_jumps);
	TR << "start: "<<(*scorefxn_)(pose) << std::endl;
	protocols::docking::DockMCMProtocolOP docker = new protocols::docking::DockMCMProtocol(1, docking_scorefxn_high_, scorefxn_);
	
	docker->set_first_cycle(first_cycle);
	docker->set_second_cycle(second_cycle);
	
	docker->apply(pose);
	TR << "end: "<<(*scorefxn_)(pose) << std::endl;
	pose.fold_tree(original_ft);
}

void
AntibodyDesignModeler::dock_low_res(core::pose::Pose& pose, std::string dock_chains, bool pack_interface /* false */) const {

	vector1< int > movable_jumps(1, 1);
	core::kinematics::FoldTree original_ft = pose.fold_tree();
	protocols::docking::setup_foldtree(pose, dock_chains, movable_jumps);
	
	TR << "start: "<<(*scorefxn_)(pose) << std::endl;
	protocols::docking::DockingLowResOP docker = new protocols::docking::DockingLowRes();
	protocols::simple_moves::SwitchResidueTypeSetMover cen_switch = protocols::simple_moves::SwitchResidueTypeSetMover("centroid");
	protocols::simple_moves::ReturnSidechainMover return_sc = protocols::simple_moves::ReturnSidechainMover(pose);
	
	cen_switch.apply(pose);
	docker->apply(pose);
	return_sc.apply(pose);
	
	pose.fold_tree(original_ft);
	
	if (pack_interface){
		repack_antigen_ab_interface(pose);
	}
	TR << "end: "<<(*scorefxn_)(pose) << std::endl;
}


}//design
}//antibody
}//protocols
