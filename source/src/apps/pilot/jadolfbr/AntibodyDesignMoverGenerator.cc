// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody_design/AntibodyDesignMoverGenerator.cc
/// @brief Handles modeling of the antibody.  Before and after design
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

//Antibody Headers
#include <protocols/antibody/design/AntibodyDesignMoverGenerator.hh>
#include <protocols/antibody/design/AntibodySeqDesignTFCreator.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>
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
#include <core/optimization/MinimizerOptions.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/kinematics/MoveMap.hh>

//Protocol Headers
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/ChangeAndResetFoldTreeMover.hh>
#include <protocols/moves/ChangeFoldTreeMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/CentroidRelax.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/docking/DockingLowRes.hh>
#include <protocols/docking/DockMCMProtocol.hh>
#include <protocols/docking/util.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <protocols/toolbox/task_operations/PreventChainFromRepackingOperation.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/loops/util.hh>
#include <protocols/loops/loops_main.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
//#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/file/FileName.hh>
#include <sstream>

static thread_local basic::Tracer TR("antibody.design.AntibodyDesignMoverGenerator");
namespace protocols {
namespace antibody {
namespace design {


using namespace core::kinematics;
using namespace core::pose;
using namespace core::scoring;
using namespace core::pack::task;
using namespace core::pack::task::operation;
using namespace protocols::simple_moves;
using namespace protocols::toolbox::task_operations;
using namespace protocols::antibody::clusters;
using namespace protocols::relax;
using namespace protocols::docking;
using core::Size;
using std::string;

AntibodyDesignMoverGenerator::AntibodyDesignMoverGenerator():
	protocols::moves::MoverApplyingMover("AntibodyDesignMoverGenerator"),
	ab_info_( /* NULL */ ),
	mover_( /* NULL */ ),
	set_tf_( /* NULL */ )
{
	ft_mover_ = protocols::moves::ChangeFoldTreeMoverOP( new protocols::moves::ChangeFoldTreeMover() );
	set_defaults();
	read_command_line_options();
	setup_scorefxns();
	setup_task_operations();

}

AntibodyDesignMoverGenerator::AntibodyDesignMoverGenerator(AntibodyInfoCOP ab_info):
	protocols::moves::MoverApplyingMover("AntibodyDesignMoverGenerator"),
	ab_info_( ab_info ),
	mover_( /* NULL */ ),
	set_tf_( /* NULL */ )
{
	ft_mover_ = protocols::moves::ChangeFoldTreeMoverOP( new protocols::moves::ChangeFoldTreeMover() );
	set_defaults();
	read_command_line_options();
	setup_scorefxns();
	setup_task_operations();

}

simple_moves::ChangeAndResetFoldTreeMoverOP
AntibodyDesignMoverGenerator::generate_protocol_mover() const {
	simple_moves::ChangeAndResetFoldTreeMoverOP protocol_mover( new simple_moves::ChangeAndResetFoldTreeMover(mover_, ft_mover_, scorefxn_) );
	return protocol_mover;
}

void
AntibodyDesignMoverGenerator::apply( core::pose::Pose & pose ){

	if ( ! mover_ ){
		utility_exit_with_message("Cannot apply mover - mover not set!!" );
	}

	core::kinematics::FoldTree original_ft = pose.fold_tree();
	//TR << "original: " << original_ft << std::endl;
	if ( use_ft_mover_ && !cartmin_){
		//Add any cutpoint variants needed. Note that the core::pose safely_add_variants does not work here (docking problem)
		ft_mover_->apply(pose);
		protocols::loops::add_cutpoint_variants(pose);
		//TR << "new: "<<pose.fold_tree() << " for "<< mover_->get_name()<< std::endl;
	}

	core::Real start_score = scorefxn_->score(pose);
	mover_->apply(pose);
	TR << "applied " << mover_->get_name() << ": start score: " << start_score << std::endl;

	TR << "applied " << mover_->get_name() << ": end   score: " << scorefxn_->score(pose) << std::endl;

	pose.fold_tree(original_ft);
	if (! cartmin_){
		protocols::loops::remove_cutpoint_variants(pose); //Remove any cutpoint variants.
	}
	if (start_coord_csts_){
		core::scoring::constraints::remove_constraints_of_type(pose, "CoordinateConstraint");
	}

	//Regenerate the FT mover.  Should be a better way than this...
	if ( !ft_mover_) {
		ft_mover_ = protocols::moves::ChangeFoldTreeMoverOP( new protocols::moves::ChangeFoldTreeMover() );
	}
}


void
AntibodyDesignMoverGenerator::set_defaults(){

	model_cdrs_.clear();
	model_cdrs_.resize(CDRNameEnum_total, true);

	overhang_ = 0;
	high_res_dock_outer_cycles_ = 3;
	high_res_dock_inner_cycles_ = 10;
	min_interface_ = false;
	min_sc_ = true;
	include_neighbor_sc_ = true;
	cartmin_ =  false;
	dualspace_ = false;
	start_coord_csts_ = false;
	set_as_mover_ = true;
	use_ft_mover_ = true;

}

AntibodyDesignMoverGenerator::~AntibodyDesignMoverGenerator(){}

void
AntibodyDesignMoverGenerator::read_command_line_options(){
	interface_detection_dis(basic::options::option [basic::options::OptionKeys::antibody::design::interface_dis]());
	neighbor_detection_dis(basic::options::option [basic::options::OptionKeys::antibody::design::neighbor_dis]());
	//set_constrain_dock_to_ab_chain(basic::options::option [basic::options::OptionKeys::antibody::design::ab_dock_chains]());
}

void
AntibodyDesignMoverGenerator::setup_scorefxns(){

	scorefxn_ = get_score_function();
	scorefxn_->apply_patch_from_file("antibody_design");
	docking_scorefxn_high_ = ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
	docking_scorefxn_high_->apply_patch_from_file("antibody_design");

	docking_scorefxn_low_ = ScoreFunctionFactory::create_score_function("interchain_cen");
	docking_scorefxn_low_->apply_patch_from_file("antibody_design");
}

void
AntibodyDesignMoverGenerator::setup_task_operations(){

	cmd_line_operation_ = core::pack::task::operation::InitializeFromCommandlineOP( new core::pack::task::operation::InitializeFromCommandline() );
	restrict_design_operation_ = core::pack::task::operation::RestrictToRepackingOP( new core::pack::task::operation::RestrictToRepacking() );
	loops_operation_ = protocols::toolbox::task_operations::RestrictToLoopsAndNeighborsOP( new RestrictToLoopsAndNeighbors() );

	tf_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory() );
	interface_tf_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory() );
	interface_tf_->push_back(cmd_line_operation_);
	interface_tf_->push_back( restrict_design_operation_); // For now.
	interface_tf_->push_back(TaskOperationCOP( new RestrictToInterface(1, interface_dis_/*Distance*/) ));
}

core::pack::task::TaskFactoryOP
AntibodyDesignMoverGenerator::get_dock_tf() {
	core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory() );
	tf->push_back(cmd_line_operation_);
	tf->push_back( restrict_design_operation_);
	tf->push_back(TaskOperationCOP( new RestrictToInterface(1, interface_dis_) ));
	return tf;
}

void
AntibodyDesignMoverGenerator::set_as_mover(bool setting) {
	set_as_mover_ = setting;
}

void
AntibodyDesignMoverGenerator::set_cdr_range(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, bool setting) {
	for (core::SSize i=cdr_start; i<=cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		set_cdr(cdr, setting);
	}
}

void
AntibodyDesignMoverGenerator::set_cdrs(utility::vector1<bool> cdrs){
	model_cdrs_ = cdrs;
	for (core::Size i = 1; i <= core::Size( CDRNameEnum_total ); ++i ){
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
		if (! ab_info_->has_CDR( cdr ) ){
			model_cdrs_[ i ] = false;
		}
	}
}


void
AntibodyDesignMoverGenerator::set_cdr(CDRNameEnum const cdr, bool setting){
	if ( ab_info_->has_CDR( cdr ) ){
		model_cdrs_[ cdr ] = setting;
	}
	else if ( setting == true ){
		TR <<"Could not set model CDR - " << ab_info_->get_CDR_name( cdr ) << " not present"<<std::endl;
		model_cdrs_[ cdr ] = false;
	}

}

void
AntibodyDesignMoverGenerator::set_cdr_only( CDRNameEnum cdr, bool setting ){
	if (setting == true){
		set_cdr_range( CDRNameEnum_start, CDRNameEnum_total, false );
		set_cdr( cdr, true );
	}
	else{
		set_cdr_range( CDRNameEnum_start, CDRNameEnum_total, true );
		set_cdr(cdr, false );
	}
}

void
AntibodyDesignMoverGenerator::stem_size(core::Size overhang) {
	overhang_ = overhang;
}

void
AntibodyDesignMoverGenerator::ab_dock_chains(std::string ab_dock_chains) {
	ab_dock_chains_ = ab_dock_chains;

}

void
AntibodyDesignMoverGenerator::set_scorefunction(ScoreFunctionCOP scorefxn){
	scorefxn_ = scorefxn->clone();
}

void
AntibodyDesignMoverGenerator::set_dock_low_res_scorefunction(ScoreFunctionCOP scorefxn){
	docking_scorefxn_low_ = scorefxn->clone();
}

void
AntibodyDesignMoverGenerator::set_dock_high_res_scorefunction(ScoreFunctionCOP scorefxn){
	docking_scorefxn_high_ = scorefxn->clone();
}

void
AntibodyDesignMoverGenerator::set_high_res_dock_inner_cycles(core::Size second_cycle /* 45 */){
	high_res_dock_inner_cycles_ = second_cycle;
}

void
AntibodyDesignMoverGenerator::set_high_res_dock_outer_cycles(core::Size first_cycle /* 4 */) {
	high_res_dock_outer_cycles_ = first_cycle;
}

void
AntibodyDesignMoverGenerator::interface_detection_dis(core::Real interface_distance){
	interface_dis_ = interface_distance;
}

void
AntibodyDesignMoverGenerator::neighbor_detection_dis(core::Real neighbor_distance){
	neighbor_dis_ = neighbor_distance;
}

void
AntibodyDesignMoverGenerator::set_task_factory(TaskFactoryOP tf){
	set_tf_ = tf;
}

void
AntibodyDesignMoverGenerator::set_mover(protocols::moves::MoverOP mover){
	mover_ = mover;
}

protocols::moves::MoverOP
AntibodyDesignMoverGenerator::mover() const {
	return mover_->clone();
}

void
AntibodyDesignMoverGenerator::set_ft_mover(protocols::moves::ChangeFoldTreeMoverOP ft_mover){
	ft_mover_ = ft_mover;
}

moves::ChangeFoldTreeMoverCOP
AntibodyDesignMoverGenerator::ft_mover() const {
	return ft_mover_;
}


void
AntibodyDesignMoverGenerator::setup_general_min_foldtree(core::pose::Pose const & pose, moves::ChangeFoldTreeMover & ft_mover){

	core::Size total_cdrs_min = 0;
	for (core::Size i = 1; i <= 6; ++i){
		if (model_cdrs_[i]) total_cdrs_min+=1;
	}

	if (min_interface_ && (total_cdrs_min > 0)){
		utility_exit_with_message("Interface and loop minimization not yet supported.");
	}
	else if (min_interface_){
		TR << "Setting RB interface FoldTree" << std::endl;
		core::kinematics::FoldTree ft;
		vector1< int > movable_jumps(1, 1);
		protocols::docking::setup_foldtree(pose, protocols::antibody::design::get_dock_chains_from_ab_dock_chains(ab_info_, ab_dock_chains_), movable_jumps, ft);
		ft_mover.set_foldtree(ft);
	}
	else {
		TR << "Setting Loop FoldTree" << std::endl;
		core::kinematics::FoldTree ft = core::kinematics::FoldTree(pose.total_residue());
		protocols::loops::LoopsOP loops = protocols::antibody::get_cdr_loops(ab_info_, pose, model_cdrs_, overhang_);
		//TR << *loops << std::endl;
		protocols::loops::fold_tree_from_loops(pose, *loops, ft);
		ft_mover.set_foldtree(ft);

	}
	TR << "Foldtree set" <<std::endl;

}

void
AntibodyDesignMoverGenerator::setup_dock_foldtree(core::pose::Pose const & pose, moves::ChangeFoldTreeMover& ft_mover){
	std::string dock_chains = get_dock_chains_from_ab_dock_chains(ab_info_, ab_dock_chains_);
	TR << "Settting up Docking FoldTree " << dock_chains << std::endl;
	vector1< int > movable_jumps(1, 1);

	core::kinematics::FoldTree ft;
	protocols::docking::setup_foldtree(pose, dock_chains, movable_jumps, ft);
	ft_mover.set_foldtree(ft);

}

////////////////////////////////////////////////////////// Modeling Functions /////////////////////////////////////////////


PackRotamersMoverOP
AntibodyDesignMoverGenerator::generate_repack_cdrs(Pose const & pose) {
	using namespace protocols::toolbox::task_operations;

	PackRotamersMoverOP packer( new PackRotamersMover(scorefxn_) );
	setup_repack_cdrs(pose, packer);
	return packer;

}

void
AntibodyDesignMoverGenerator::setup_repack_cdrs(Pose const & pose, PackRotamersMoverOP packer) {
	if (set_tf_){
		packer->task_factory(set_tf_);
	}
	else {
		protocols::loops::LoopsOP cdr_loops = protocols::antibody::get_cdr_loops(ab_info_, pose, model_cdrs_, overhang_);
		loops_operation_->set_design_loop(false);
		loops_operation_->set_include_neighbors(include_neighbor_sc_);
		loops_operation_->set_cutoff_distance(neighbor_dis_);
		loops_operation_->set_loops(cdr_loops);

		tf_->clear();
		tf_->push_back(cmd_line_operation_ );
		tf_->push_back( restrict_design_operation_ );
		tf_->push_back(loops_operation_);
		packer->task_factory(tf_);
	}

	if (set_as_mover_){
		use_ft_mover_ = false;
		mover_ = packer;
	}
}

PackRotamersMoverOP
AntibodyDesignMoverGenerator::generate_repack_antigen_ab_interface( Pose const & pose ){
	vector1< char > antigen_chains = ab_info_->get_antigen_chains();
	if (antigen_chains.size() == 0){
		TR <<" Antigen not present to repack interface!" << std::endl;
		return NULL;
	}

	PackRotamersMoverOP packer( new PackRotamersMover() );
	setup_repack_antigen_ab_interface(pose, packer);
	return packer;
}

void
AntibodyDesignMoverGenerator::setup_repack_antigen_ab_interface(Pose const& pose, simple_moves::PackRotamersMoverOP packer) {
	vector1< char > antigen_chains = ab_info_->get_antigen_chains();
	if (antigen_chains.size() == 0){
		TR <<" Antigen not present to repack interface!" << std::endl;
		return;
	}

	packer->score_function(scorefxn_);
	if (set_tf_){
		packer->task_factory(set_tf_);
	}
	else {
		packer->task_factory(interface_tf_);
	}

	if (set_as_mover_){
		std::string current_dock_chains = ab_dock_chains_;
		ab_dock_chains_ = "LH_A";
		this->setup_dock_foldtree( pose, *ft_mover_ );
		use_ft_mover_ = true;
		mover_ = packer;
		ab_dock_chains_ = current_dock_chains;
	}
}

MinMoverOP
AntibodyDesignMoverGenerator::generate_minimizer(Pose & pose) {

	//Tolerance of .001 is slower by ~25% or less.  However, it can also decrease energies significantly compared to .01.
	//(In a test set of one protein minimized by bb - write a script to test this )
	MinMoverOP min_mover( new MinMover() );
	setup_minimizer( pose, min_mover);
	return min_mover;

}

void
AntibodyDesignMoverGenerator::setup_minimizer(Pose& pose, MinMoverOP min_mover){

	ScoreFunctionOP local_scorefxn;
	if (cartmin_){
		local_scorefxn = scorefxn_->clone();
		local_scorefxn->set_weight_if_zero(cart_bonded, .5);
		local_scorefxn->set_weight(pro_close, 0);
	} else {
		local_scorefxn =  scorefxn_->clone();
	}

	local_scorefxn->set_weight(chainbreak, 100.00);

	MoveMapOP mm;
	if (! min_interface_ ){
		mm = get_cdrs_movemap_with_overhang(pose, true /*BB*/, min_sc_/*SC*/, include_neighbor_sc_, false /*neighbor_bb*/);
	}

	min_mover->set_movemap(mm);
	min_mover->score_function(local_scorefxn);
	min_mover->min_type("dfpmin_armijo_nonmonotone");
	min_mover->tolerance(0.001);

	if (cartmin_){
		min_mover->min_type("lbfgs_armijo_nonmonotone");
		min_mover->cartesian(true);
		min_mover->min_options()->max_iter(200); //Great suggestion from Patrick Conway
	}

	if (set_as_mover_){

		if ( cartmin_ ){
			use_ft_mover_ = false;
		}
		else{
			this->setup_general_min_foldtree( pose, *ft_mover_ );
			use_ft_mover_ = true;
		}
		mover_ = min_mover;
	}
}

FastRelaxOP
AntibodyDesignMoverGenerator::generate_relax( Pose & pose ){


	protocols::relax::FastRelaxOP rel( new protocols::relax::FastRelax() );//Stack construction failed to compile.
	setup_relax(pose, rel);
	return rel;
}

void
AntibodyDesignMoverGenerator::setup_relax(Pose & pose, relax::FastRelaxOP rel){

	ScoreFunctionOP local_scorefxn = scorefxn_->clone();
	local_scorefxn->set_weight(chainbreak, 100.00);

	protocols::loops::LoopsOP cdr_loops = get_cdr_loops(ab_info_, pose, model_cdrs_, overhang_);

	if ( set_tf_ ){
		rel->set_task_factory(set_tf_);
	}
	else {
		loops_operation_->set_design_loop(false);
		loops_operation_->set_include_neighbors(include_neighbor_sc_);
		loops_operation_->set_loops(cdr_loops);
		loops_operation_->set_cutoff_distance(neighbor_dis_);

		tf_->clear();
		tf_->push_back(cmd_line_operation_ );
		tf_->push_back( restrict_design_operation_ );
		tf_->push_back(loops_operation_);
		rel->set_task_factory(tf_);
	}

	MoveMapOP mm = get_cdrs_movemap_with_overhang(pose, true /*BB*/, min_sc_/*SC*/, include_neighbor_sc_,  false /*neighbor_bb*/);

	if (dualspace_) {
		local_scorefxn->set_weight_if_zero( cart_bonded, .5 );
		local_scorefxn->set_weight( pro_close, 0 );
		rel->dualspace( true );
		rel->max_iter( 200 );
		//rel->minimize_bond_angles(true); Too Long - should be option
	}
	rel->set_scorefxn( local_scorefxn );
	rel->set_movemap( mm );

	if ( start_coord_csts_ ){
		rel->constrain_relax_to_start_coords( true );
	}

	if ( set_as_mover_ ){
		cartmin_ = false; //This mover is way to complicated for its own good.  Seriously.  Should have went with separate setup movers.
		this->setup_general_min_foldtree( pose, *ft_mover_ );
		use_ft_mover_ = true;
		mover_ = rel;

	}

}

DockingLowResOP
AntibodyDesignMoverGenerator::generate_dock_low_res( Pose const & pose ) {

	protocols::docking::DockingLowResOP docker( new protocols::docking::DockingLowRes() );
	setup_dock_low_res(pose, docker);
	return docker;
}

void
AntibodyDesignMoverGenerator::setup_dock_low_res(const Pose& pose, DockingLowResOP docker) {
	docker->set_scorefxn(docking_scorefxn_low_);
	protocols::simple_moves::SwitchResidueTypeSetMoverOP cen_switch( new protocols::simple_moves::SwitchResidueTypeSetMover("centroid") );
	protocols::simple_moves::ReturnSidechainMoverOP return_sc( new protocols::simple_moves::ReturnSidechainMover(pose) );

	if ( set_as_mover_ ){
		protocols::moves::SequenceMoverOP seq_mover( new moves::SequenceMover() );
		seq_mover->add_mover(cen_switch);
		seq_mover->add_mover(docker);
		seq_mover->add_mover(return_sc);
		mover_ = seq_mover;
		this->setup_dock_foldtree( pose, *ft_mover_ );
		use_ft_mover_ = true;
	}
}

protocols::docking::DockMCMProtocolOP
AntibodyDesignMoverGenerator::generate_dock_high_res( Pose const & pose ) {

	protocols::docking::DockMCMProtocolOP docker( new protocols::docking::DockMCMProtocol(1, docking_scorefxn_high_, scorefxn_) );
	setup_dock_high_res(pose, docker);

	return docker;
}

void
AntibodyDesignMoverGenerator::setup_dock_high_res(const Pose& pose, DockMCMProtocolOP docker){

	docker->set_first_cycle(high_res_dock_outer_cycles_);
	docker->set_second_cycle(high_res_dock_inner_cycles_);
	docker->set_scorefxn(docking_scorefxn_high_);
	docker->set_scorefxn_pack(scorefxn_);

	if (set_tf_){
		docker->set_task_factory(set_tf_);
		docker->set_ignore_default_task(true);

		//Debuging
		//core::pack::task::PackerTaskOP task = set_tf_->create_task_and_apply_taskoperations(pose);
		//task->show();
	}
	else{
		docker->set_task_factory(this->get_dock_tf());
		docker->set_ignore_default_task(true);
	}

	if ( set_as_mover_ ){
		mover_ = docker;
		this->setup_dock_foldtree( pose, *ft_mover_ );
		use_ft_mover_ = true;
	}

}

MoveMapOP
AntibodyDesignMoverGenerator::get_cdrs_movemap_with_overhang(Pose  & pose, bool min_bb, bool min_sc, bool include_neighbor_sc, bool include_neighbor_bb) const {

	pose.update_residue_neighbors();
	MoveMapOP mm( new MoveMap() );
	for (core::Size i = 1; i <= core::Size( ab_info_->get_total_num_CDRs() ); ++i){
		if (model_cdrs_[i]){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			protocols::loops::Loop cdr_loop = ab_info_->get_CDR_loop(cdr, pose, overhang_);

			vector1<bool> all_included_res(pose.total_residue(), false);
			select_loop_residues( pose, cdr_loop, (include_neighbor_sc || include_neighbor_bb) ,all_included_res , neighbor_dis_);
			for (core::Size x = 1; x <= pose.total_residue(); ++x){
				if (x >= cdr_loop.start() && x <= cdr_loop.stop()){
					if (min_bb) mm->set_bb(x, true);
					if (min_sc) mm->set_chi(x, true);

				}
				//Not in the loop, then it's a neighbor
				else if  (all_included_res[x]) {
					if (include_neighbor_sc) mm->set_chi(x, true);
					if (include_neighbor_bb) mm->set_bb(x, true);
				}
			}
		}
	}

	return mm;
}

MoveMapOP
AntibodyDesignMoverGenerator::get_movemap_from_task(core::pose::Pose const & pose, core::pack::task::PackerTaskCOP task) const {
	MoveMapOP mm( new MoveMap() );
	for (core::Size i = 1; i<=pose.total_residue(); i++){
		if (task->pack_residue(i)){
			mm->set_chi(i, true);
		}
	}
	return mm;
}

void
AntibodyDesignMoverGenerator::set_min_interface( bool min_interface ){ min_interface_ = min_interface; }

void
AntibodyDesignMoverGenerator::set_min_sc( bool min_sc ) { min_sc_ = min_sc; }

void
AntibodyDesignMoverGenerator::set_include_neighbor_sc( bool include_neighbor_sc ) { include_neighbor_sc_ = include_neighbor_sc; }

void
AntibodyDesignMoverGenerator::set_cartmin( bool cartmin ) { cartmin_ = cartmin; }

void
AntibodyDesignMoverGenerator::set_dualspace( bool dualspace ) { dualspace_ = dualspace; }

void
AntibodyDesignMoverGenerator::set_start_coord_csts( bool coord_csts ) { start_coord_csts_ = coord_csts; }

}//design
}//antibody
}//protocols

