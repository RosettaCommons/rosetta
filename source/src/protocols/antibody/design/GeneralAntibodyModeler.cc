// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody_design/GeneralAntibodyModeler.cc
/// @brief Handles modeling of the antibody.  Before and after design
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

//Antibody Headers
#include <protocols/antibody/design/GeneralAntibodyModeler.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>
#include <protocols/antibody/AntibodyInfo.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/chains_util.hh>
#include <core/scoring/ScoreType.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/util.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/kinematics/MoveMap.fwd.hh>

//Protocol Headers
#include <protocols/backrub/BackrubProtocol.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/docking/DockingLowRes.hh>
#include <protocols/docking/DockMCMProtocol.hh>
#include <protocols/docking/util.hh>
#include <protocols/simple_task_operations/RestrictToInterface.hh>
#include <protocols/task_operations/PreventChainFromRepackingOperation.hh>
#include <protocols/simple_task_operations/RestrictToLoopsAndNeighbors.hh>
#include <protocols/loops/util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loops.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

//#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility>
#include <sstream>

#include <protocols/simple_moves/ReturnSidechainMover.hh> // AUTO IWYU For ReturnSidechainMover

static basic::Tracer TR( "antibody.design.GeneralAntibodyModeler" );
namespace protocols {
namespace antibody {
namespace design {


using namespace core::kinematics;
using namespace core::pose;
using namespace core::scoring;
using namespace core::pack::task;
using namespace core::pack::task::operation;
using namespace protocols::simple_moves;
using namespace protocols::minimization_packing;
using namespace protocols::task_operations;
using namespace protocols::simple_task_operations;
using namespace protocols::antibody::clusters;
using core::Size;
using std::string;
using utility::vector1;

GeneralAntibodyModeler::GeneralAntibodyModeler(AntibodyInfoOP ab_info) :
	utility::VirtualBase(),
	ab_info_(std::move(ab_info))
{


	set_defaults();
	read_command_line_options();
	setup_scorefxns();
	setup_task_operations();

	core::pose::DockingPartners dock_chains;
	dock_chains.partner1.push_back("A");
	dock_chains.partner2 = ab_info_->get_antibody_chains();
	ab_dock_chains(dock_chains);

}

void
GeneralAntibodyModeler::set_defaults(){

	model_cdrs_.clear();
	model_cdrs_.resize(CDRNameEnum_proto_total, true);
	model_cdrs_[ l4 ] = false;
	model_cdrs_[ h4 ] = false;

	overhangs_.clear();
	overhangs_.resize(CDRNameEnum_proto_total, 0);

}

GeneralAntibodyModeler::~GeneralAntibodyModeler()= default;

GeneralAntibodyModeler::GeneralAntibodyModeler( GeneralAntibodyModeler const & src ):
	utility::VirtualBase(),
	interface_dis_( src.interface_dis_),
	neighbor_dis_( src.neighbor_dis_),
	ab_dock_chains_(src.ab_dock_chains_),
	model_cdrs_( src.model_cdrs_ ),
	overhangs_( src.overhangs_ ),
	tf_set_( src.tf_set_),
	atom_pair_weight_( src.atom_pair_weight_ )

{
	using namespace protocols::task_operations;
	using namespace core::pack::task::operation;
	if ( src.ab_info_ ) ab_info_ = src.ab_info_->clone();
	if ( src.scorefxn_ ) scorefxn_ = src.scorefxn_->clone();
	if ( src.min_scorefxn_ ) min_scorefxn_ = src.min_scorefxn_->clone();
	if ( src.docking_scorefxn_low_ ) docking_scorefxn_low_ = src.docking_scorefxn_low_->clone();
	if ( src.docking_scorefxn_high_ ) docking_scorefxn_high_ = src.docking_scorefxn_high_->clone();

	if ( src.interface_tf_ ) interface_tf_ = src.interface_tf_->clone();
	if ( src.tf_ ) tf_ = src.tf_->clone();

	if ( src.loops_operation_ ) {
		loops_operation_ = utility::pointer::make_shared< RestrictToLoopsAndNeighbors >( *src.loops_operation_ );
	}
	if ( src.cmd_line_operation_ ) {
		cmd_line_operation_ = utility::pointer::make_shared< InitializeFromCommandline >( *src.cmd_line_operation_ );
	}

	if ( src.restrict_design_operation_ ) {
		restrict_design_operation_ = utility::pointer::make_shared< RestrictToRepacking >( *src.restrict_design_operation_ );
	}



}

GeneralAntibodyModelerOP
GeneralAntibodyModeler::clone() const {
	return utility::pointer::make_shared< GeneralAntibodyModeler >( *this );
}
void
GeneralAntibodyModeler::read_command_line_options(){

	interface_detection_dis(basic::options::option [basic::options::OptionKeys::antibody::design::interface_dis]());
	neighbor_detection_dis(basic::options::option [basic::options::OptionKeys::antibody::design::neighbor_dis]());
	atom_pair_weight_ = basic::options::option [ basic::options::OptionKeys::antibody::design::atom_pair_cst_weight]();
}

void
GeneralAntibodyModeler::setup_scorefxns(){

	using namespace basic::options;

	scorefxn_ = design::get_ab_design_global_scorefxn();
	min_scorefxn_ = design::get_ab_design_min_scorefxn();
	docking_scorefxn_high_ = design::get_ab_design_dock_high_scorefxn();
	docking_scorefxn_low_ = design::get_ab_design_dock_low_scorefxn();
}

void
GeneralAntibodyModeler::setup_task_operations(){

	tf_set_ = false;
	cmd_line_operation_ = utility::pointer::make_shared< core::pack::task::operation::InitializeFromCommandline >();
	restrict_design_operation_ = utility::pointer::make_shared< core::pack::task::operation::RestrictToRepacking >();
	loops_operation_ = utility::pointer::make_shared< RestrictToLoopsAndNeighbors >();

	tf_ = utility::pointer::make_shared< core::pack::task::TaskFactory >();
	interface_tf_ = utility::pointer::make_shared< core::pack::task::TaskFactory >();
	interface_tf_->push_back(cmd_line_operation_);
	interface_tf_->push_back( restrict_design_operation_);
	interface_tf_->push_back(utility::pointer::make_shared< RestrictToInterface >(1, interface_dis_/*Distance*/));

	interface_tf_->push_back( utility::pointer::make_shared< ExtraRotamers >( 0 /*all*/, 1 /*ex1*/, 1 /*level*/ ) );
	interface_tf_->push_back( utility::pointer::make_shared< ExtraRotamers >( 0 /*all*/, 2 /*ex2*/, 1 /*level*/ ) );
	interface_tf_->push_back(utility::pointer::make_shared< IncludeCurrent >());
}

void
GeneralAntibodyModeler::set_cdr_range(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, bool setting) {
	for ( core::SSize i=cdr_start; i<=cdr_end; ++i ) {
		auto cdr = static_cast<CDRNameEnum>(i);
		set_cdr(cdr, setting);
	}
}

void
GeneralAntibodyModeler::set_cdrs(utility::vector1<bool> cdrs){
	model_cdrs_ = cdrs;
	if ( model_cdrs_.size() == 6 ) {
		model_cdrs_.push_back(false);
		model_cdrs_.push_back(false);
	}

}

void
GeneralAntibodyModeler::set_cdr_only(CDRNameEnum cdr, bool setting){
	if ( setting==true ) {
		set_cdr_range(CDRNameEnum_start, CDRNameEnum_proto_total, false);
		set_cdr(cdr, true);
	} else {
		set_cdr_range(CDRNameEnum_start, CDRNameEnum_proto_total, true);
		set_cdr(cdr, false);
	}
}

void
GeneralAntibodyModeler::set_cdr(CDRNameEnum const cdr, bool setting){
	model_cdrs_[cdr]=setting;
}

void
GeneralAntibodyModeler::set_overhang( const core::Size overhang ){

	overhangs_.clear();
	overhangs_.resize(CDRNameEnum_proto_total, overhang);

}

void
GeneralAntibodyModeler::cdr_overhang(const CDRNameEnum cdr, const core::Size overhang) {
	overhangs_[cdr] = overhang;
}

//void
//GeneralAntibodyModeler::ab_dock_chains(std::string const & ab_dock_chains) {
// ab_dock_chains_ = core::pose::docking_partners_from_string(ab_dock_chains);
//}

void
GeneralAntibodyModeler::ab_dock_chains(core::pose::DockingPartners const & ab_dock_chains) {
	ab_dock_chains_ = ab_dock_chains;
}

void
GeneralAntibodyModeler::set_scorefunction(ScoreFunctionCOP scorefxn){
	scorefxn_ = scorefxn->clone();
}

void
GeneralAntibodyModeler::set_scorefunction_min(ScoreFunctionCOP scorefxn){
	min_scorefxn_ = scorefxn->clone();
}

void
GeneralAntibodyModeler::set_dock_low_res_scorefunction(ScoreFunctionCOP scorefxn){
	docking_scorefxn_low_ = scorefxn->clone();
}

void
GeneralAntibodyModeler::set_dock_high_res_scorefunction(ScoreFunctionCOP scorefxn){
	docking_scorefxn_high_ = scorefxn->clone();
}

void
GeneralAntibodyModeler::set_task_factory(core::pack::task::TaskFactoryOP tf) {
	tf_ = tf;
	tf_set_ = true;
}

void
GeneralAntibodyModeler::reset_task_factory() {
	setup_task_operations();
}

void
GeneralAntibodyModeler::interface_detection_dis(core::Real interface_distance){
	interface_dis_ = interface_distance;
}

void
GeneralAntibodyModeler::neighbor_detection_dis(core::Real neighbor_distance){
	neighbor_dis_ = neighbor_distance;
}

protocols::loops::LoopsOP
GeneralAntibodyModeler::get_cdr_loops(Pose const & pose) const {

	protocols::loops::LoopsOP cdr_loops( new protocols::loops::Loops );

	for ( auto const & cdr : ab_info_->get_all_cdrs_present( true /* include_cdr4 */) ) {
		if ( model_cdrs_[ cdr ] ) {
			core::Size start = ab_info_->get_CDR_start(cdr, pose);
			core::Size stop =  ab_info_->get_CDR_end(cdr, pose);
			core::Size cutpoint = (stop-start)/2+start;
			protocols::loops::Loop cdr_loop = protocols::loops::Loop(start, stop, cutpoint);
			cdr_loops->add_loop(cdr_loop);
		}
	}
	return cdr_loops;
}

protocols::loops::LoopsOP
GeneralAntibodyModeler::get_cdr_loops_with_overhang(Pose const & pose) const {

	protocols::loops::LoopsOP cdr_loops( new protocols::loops::Loops );

	for ( auto const & cdr : ab_info_->get_all_cdrs_present( true /* include_cdr4 */) ) {
		if ( model_cdrs_[cdr] ) {
			protocols::loops::Loop cdr_loop = get_cdr_loop_with_overhang(pose, cdr);
			cdr_loops->add_loop(cdr_loop);
		}
	}
	return cdr_loops;
}

protocols::loops::Loop
GeneralAntibodyModeler::get_cdr_loop_with_overhang(Pose const & pose, CDRNameEnum cdr) const {

	core::Size start = ab_info_->get_CDR_start(cdr, pose) - overhangs_[cdr];
	core::Size stop =  ab_info_->get_CDR_end(cdr, pose) + overhangs_[cdr];
	core::Size cutpoint = (stop-start)/2+start;
	protocols::loops::Loop cdr_loop = protocols::loops::Loop(start, stop, cutpoint);
	return cdr_loop;

}

void
GeneralAntibodyModeler::apply_A_LH_foldtree(core::pose::Pose & pose) const {
	core::pose::DockingPartners dock_chains;
	dock_chains.partner1 = ab_info_->get_antigen_chains();
	dock_chains.partner2 = utility::vector1<std::string>{"L", "H"};

	vector1< int > movable_jumps(1, 1);
	protocols::docking::setup_foldtree(pose, dock_chains, movable_jumps);
}

MoveMapOP
GeneralAntibodyModeler::get_cdrs_movemap_with_overhang(Pose  & pose, bool min_bb, bool min_sc, bool include_neighbor_sc, bool include_neighbor_bb) const {


	pose.update_residue_neighbors();
	MoveMapOP mm( new MoveMap() );
	for ( auto const & cdr : ab_info_->get_all_cdrs_present( true /* include_cdr4 */) ) {
		if ( model_cdrs_[ cdr ] ) {
			protocols::loops::Loop cdr_loop = ab_info_->get_CDR_loop(cdr, pose, overhangs_[ cdr ]);

			vector1<bool> all_included_res(pose.size(), false);
			select_loop_residues( pose, cdr_loop, (include_neighbor_sc || include_neighbor_bb) ,all_included_res , neighbor_dis_);
			for ( core::Size x = 1; x <= pose.size(); ++x ) {
				if ( x >= cdr_loop.start() && x <= cdr_loop.stop() ) {
					if ( min_bb ) mm->set_bb(x, true);
					if ( min_sc ) mm->set_chi(x, true);

				} else if  ( all_included_res[x] ) {
					//Not in the loop, then it's a neighbor
					if ( include_neighbor_sc ) mm->set_chi(x, true);
					if ( include_neighbor_bb ) mm->set_bb(x, true);
				}
			}
		}
	}

	return mm;
}

MoveMapOP
GeneralAntibodyModeler::get_movemap_from_task(core::pose::Pose const & pose, core::pack::task::PackerTaskCOP task) const {
	MoveMapOP mm( new MoveMap() );
	for ( core::Size i = 1; i<=pose.size(); i++ ) {
		if ( task->pack_residue(i) ) {
			mm->set_chi(i, true);
		}
	}
	return mm;
}

core::pose::DockingPartners
GeneralAntibodyModeler::get_dock_chains() const {
	return get_dock_chains_from_ab_dock_chains(ab_info_, ab_dock_chains_);
}

////////////////////////////////////////////////////////// Modeling Functions /////////////////////////////////////////////

void
GeneralAntibodyModeler::extend_CDR(Pose& pose, CDRNameEnum cdr) const {

	core::Size start = ab_info_->get_CDR_start(cdr, pose);
	core::Size end = ab_info_->get_CDR_end(cdr, pose);

	protocols::loops::Loop cdr_loop  = protocols::loops::Loop(start+2, end-2, start+3);
	protocols::loops::Loops cdr_loops = protocols::loops::Loops();
	cdr_loops.add_loop(cdr_loop);

	core::kinematics::FoldTree original_ft = pose.fold_tree();
	core::kinematics::FoldTree loop_ft = core::kinematics::FoldTree();

	protocols::loops::fold_tree_from_loops(pose, cdr_loops, loop_ft);
	pose.fold_tree(loop_ft);
	protocols::loops::set_extended_torsions_and_idealize_loops(pose, cdr_loops);

	pose.fold_tree(original_ft);

}

void
GeneralAntibodyModeler::repack_cdrs(Pose& pose, bool include_neighbor_sc ) {
	using namespace protocols::task_operations;

	protocols::loops::LoopsOP cdr_loops = get_cdr_loops_with_overhang(pose);

	core::Real start_e = (*scorefxn_)(pose);


	if ( ! tf_set_ ) {
		loops_operation_->set_design_loop(false);
		loops_operation_->set_include_neighbors(include_neighbor_sc);
		loops_operation_->set_cutoff_distance(neighbor_dis_);
		loops_operation_->set_loops(cdr_loops);
		tf_->clear();
		tf_->push_back(cmd_line_operation_ );
		tf_->push_back( restrict_design_operation_ );
		tf_->push_back(loops_operation_);
		tf_->push_back( utility::pointer::make_shared< ExtraRotamers >( 0 /*all*/, 1 /*ex1*/, 1 /*level*/ ) );
		tf_->push_back( utility::pointer::make_shared< ExtraRotamers >( 0 /*all*/, 2 /*ex2*/, 1 /*level*/ ) );
		tf_->push_back(utility::pointer::make_shared< IncludeCurrent >());
	}

	core::pack::task::PackerTaskOP task = tf_->create_task_and_apply_taskoperations(pose);

	PackRotamersMover packer = PackRotamersMover(scorefxn_, task);
	packer.apply(pose);
	TR << "repacking cdrs start: "<<start_e << std::endl;
	TR << "repacking cdrs end: "<<(*scorefxn_)(pose) << std::endl;

}

void
GeneralAntibodyModeler::minimize_cdrs(Pose & pose,
	bool min_sc,
	bool include_neighbor_sc ,
	bool min_interface,
	bool cartmin,
	bool use_talaris_cartmin) const {

	ScoreFunctionOP local_scorefxn;
	if ( cartmin && use_talaris_cartmin ) {
		local_scorefxn = core::scoring::ScoreFunctionFactory::create_score_function("talaris2013_cart");
	} else if ( cartmin ) {
		local_scorefxn = min_scorefxn_->clone();
		local_scorefxn->set_weight(cart_bonded, .5);
		local_scorefxn->set_weight(pro_close, 0);
	} else {
		local_scorefxn =  min_scorefxn_->clone();
	}

	local_scorefxn->set_weight(chainbreak, 100.00);
	core::Real start_e = (*scorefxn_)(pose);

	core::kinematics::FoldTree original_ft = pose.fold_tree();
	if ( min_interface ) {
		TR << "Including RB interface min" <<std::endl;
		vector1< int > movable_jumps(1, 1);
		protocols::docking::setup_foldtree(pose, get_dock_chains_from_ab_dock_chains(ab_info_, ab_dock_chains_), movable_jumps);
	} else if ( ! cartmin ) {
		//Make sure Minimization of loops does not effect the rest of the chain
		protocols::loops::LoopsOP loops = get_cdr_loops_with_overhang(pose);
		core::kinematics::FoldTree ft = core::kinematics::FoldTree();
		protocols::loops::fold_tree_from_loops(pose, *loops, ft);
		pose.fold_tree(ft);
		protocols::loops::add_cutpoint_variants(pose);
	}

	MoveMapOP mm = get_cdrs_movemap_with_overhang(pose, true /*BB*/, min_sc/*SC*/, include_neighbor_sc, false /*neighbor_bb*/);

	//Tolerance of .001 is slower by ~25% or less.  However, it can also decrease energies significantly compared to .01.
	//(In a test set of one protein minimized by bb - write a script to test this )

	protocols::minimization_packing::MinMover min_mover(mm, local_scorefxn, "lbfgs_armijo_nonmonotone", 0.001, true /*use_nblist*/ );


	if ( cartmin ) {
		min_mover.min_type("lbfgs_armijo_nonmonotone");
		min_mover.cartesian(true);
		min_mover.min_options()->max_iter(200); //Great suggestion from Patrick Conway
		min_mover.apply(pose);
	} else {
		min_mover.apply(pose);
		protocols::loops::remove_cutpoint_variants(pose);
		pose.fold_tree(original_ft);
	}

	//local_scorefxn->show(TR, pose);
	TR << "CDR min start: "<<start_e << std::endl;
	TR << "CDR min end: "<<(*scorefxn_)(pose) << std::endl;
}

void
GeneralAntibodyModeler::relax_cdrs(Pose & pose,  bool include_neighbor_sc /*true*/,  bool starting_coordinate_constraints /*false*/, bool min_interface /* false */, bool dualspace /*false*/, core::Size repeats) const {

	ScoreFunctionOP local_scorefxn = min_scorefxn_->clone();


	local_scorefxn->set_weight(chainbreak, 100.00);

	core::kinematics::FoldTree original_ft = pose.fold_tree();
	protocols::loops::LoopsOP cdr_loops = get_cdr_loops_with_overhang(pose);

	if ( min_interface ) {
		TR << "Including RB interface min" <<std::endl;
		vector1< int > movable_jumps(1, 1);
		protocols::docking::setup_foldtree(pose, protocols::antibody::design::get_dock_chains_from_ab_dock_chains(ab_info_, ab_dock_chains_), movable_jumps);
	} else {
		//else if (! dualspace){
		//Make sure Minimization of loops does not effect the rest of the chain
		//Dualspace as well - since half of it is dihedral, we still need cutpoints.
		//No reason chainbreak term won't work with cart.

		core::kinematics::FoldTree ft = core::kinematics::FoldTree();
		protocols::loops::fold_tree_from_loops(pose, *cdr_loops, ft);
		pose.fold_tree(ft);
		protocols::loops::add_cutpoint_variants(pose);
	}

	if ( ! tf_set_ ) {
		loops_operation_->set_design_loop(false);
		loops_operation_->set_include_neighbors(include_neighbor_sc);
		loops_operation_->set_loops(cdr_loops);
		loops_operation_->set_cutoff_distance(neighbor_dis_);

		tf_->clear();
		tf_->push_back(cmd_line_operation_ );
		tf_->push_back( restrict_design_operation_ );
		tf_->push_back(loops_operation_);

		tf_->push_back( utility::pointer::make_shared< ExtraRotamers >( 0 /*all*/, 1 /*ex1*/, 1 /*level*/ ) );
		tf_->push_back( utility::pointer::make_shared< ExtraRotamers >( 0 /*all*/, 2 /*ex2*/, 1 /*level*/ ) );
		tf_->push_back(utility::pointer::make_shared< IncludeCurrent >());
	}

	core::Real start_e = (*scorefxn_)(pose);
	MoveMapOP mm = get_cdrs_movemap_with_overhang(pose, true /*BB*/, true/*SC*/, include_neighbor_sc,  false /*neighbor_bb*/);


	protocols::relax::FastRelaxOP rel( new protocols::relax::FastRelax( repeats ) );//Stack construction failed to compile.

	if ( dualspace ) {
		local_scorefxn->set_weight(cart_bonded, .5);
		local_scorefxn->set_weight(pro_close, 0);
		rel->dualspace(true);
		rel->max_iter(200);
		//rel->minimize_bond_angles(true); Too Long - should be option
	}

	rel->set_scorefxn(local_scorefxn);
	rel->set_movemap(mm);
	//tf_->create_task_and_apply_taskoperations(pose)->show(TR);

	rel->set_task_factory(tf_);
	if ( starting_coordinate_constraints ) {
		rel->constrain_relax_to_start_coords(true);
	}
	rel->apply(pose);

	//local_scorefxn->show_pretty(TR.Debug);

	TR << "CDR relax start: "<<start_e << std::endl;
	TR << "CDR relax end: "<<(*scorefxn_)(pose) << std::endl;
	if ( starting_coordinate_constraints ) {
		core::scoring::constraints::remove_constraints_of_type(pose, "CoordinateConstraint");
	}
	protocols::loops::remove_cutpoint_variants(pose);
	pose.fold_tree(original_ft);
}

void
GeneralAntibodyModeler::relax_interface(Pose & pose, bool min_interface_sc /* true */) const {

	//TR << "Relaxing Interface: " << dock_chains << std::endl;
	ScoreFunctionOP local_scorefxn = scorefxn_->clone();
	vector1< int > movable_jumps(1, 1);
	core::kinematics::FoldTree original_ft = pose.fold_tree();
	protocols::docking::setup_foldtree(pose, get_dock_chains_from_ab_dock_chains(ab_info_, ab_dock_chains_), movable_jumps);

	MoveMapOP mm( utility::pointer::make_shared< MoveMap >() );
	mm->set_jump(1, true); //There is more than one jump!

	core::pack::task::TaskFactoryOP tf = interface_tf_->clone();
	core::pack::task::PackerTaskOP task = interface_tf_->create_task_and_apply_taskoperations(pose);
	protocols::relax::FastRelaxOP rel( new protocols::relax::FastRelax(local_scorefxn) );

	//Add residues to movemap, so that they minimize during minimization, while repacking will happen via tf. Maybe should not have any chi packing.
	if ( min_interface_sc ) {
		mm = get_movemap_from_task(pose, task);
	}

	rel->set_movemap(mm);
	if ( !tf_set_ ) {
		rel->set_task_factory(tf);
	} else {
		rel->set_task_factory(tf_);
	}
	rel->apply(pose);
	pose.fold_tree(original_ft);
}

void
GeneralAntibodyModeler::minimize_interface(Pose& pose, bool min_interface_sc /* true */) const {


	vector1< int > movable_jumps(1, 1);
	core::kinematics::FoldTree original_ft = pose.fold_tree();
	core::pose::DockingPartners dock_chains = get_dock_chains_from_ab_dock_chains(ab_info_, ab_dock_chains_);
	TR << "Minimizing Interface: " << dock_chains << std::endl;
	protocols::docking::setup_foldtree(pose, dock_chains, movable_jumps);
	core::pack::task::PackerTaskOP task = interface_tf_->create_task_and_apply_taskoperations(pose);
	MoveMapOP mm;
	if ( min_interface_sc ) {
		mm = get_movemap_from_task(pose, task);
	}

	mm->set_jump(1, true); //There is more than one jump!

	ScoreFunctionOP local_scorefxn = scorefxn_->clone();
	local_scorefxn->set_weight_if_zero(atom_pair_constraint, atom_pair_weight_);

	protocols::minimization_packing::MinMover min_mover(mm, local_scorefxn, "lbfgs_armijo_nonmonotone", 0.01, true /*use_nblist*/ );
	min_mover.apply(pose);
	pose.fold_tree(original_ft);
}

void
GeneralAntibodyModeler::backrub_cdrs( core::pose::Pose & pose, bool min_sc, bool include_neighbor_sc /*true*/) const {

	using namespace protocols::backrub;
	using namespace basic::options;

	ScoreFunctionOP local_scorefxn  = min_scorefxn_->clone();
	//local_scorefxn->set_weight(chainbreak, 20.00);

	local_scorefxn->set_weight(dihedral_constraint, 0); //Having this not zero will probably mess up the backrub motion
	//Note that if global_atom_pair_cst_scoring is on, then this may effect backrub.  Which may be fine. I guess.

	core::Real start_e = (*scorefxn_)(pose);

	//core::kinematics::FoldTree original_ft = pose.fold_tree();


	//Make sure Minimization of loops does not effect the rest of the chain
	protocols::loops::LoopsOP loops = get_cdr_loops_with_overhang(pose);
	//core::kinematics::FoldTree ft = core::kinematics::FoldTree();
	//protocols::loops::fold_tree_from_loops(pose, *loops, ft);
	//pose.fold_tree(ft);
	//protocols::loops::add_cutpoint_variants(pose);

	MoveMapOP mm = get_cdrs_movemap_with_overhang(pose, true /*BB*/, min_sc/*SC*/, include_neighbor_sc, false /*neighbor_bb*/);

	//Additional pivots should be at least 3 residues long.  We need to pass these as additional loops to RestrictToLoopsAndNeighbors so we pack these.
	// As well as add them to the movemap
	if ( option[ OptionKeys::antibody::design::add_backrub_pivots].user() ) {
		utility::vector1<std::string> res_strings = option [ OptionKeys::antibody::design::add_backrub_pivots]();
		utility::vector1<bool> add_residues = get_resnums_from_strings_with_ranges(pose, res_strings);
		add_loops_from_bool_vector(*loops, add_residues);
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			if ( add_residues[ i ] &&  i >= 1 ) {
				mm->set_bb( i , true);
				mm->set_chi(i, true);
			}

		}
	}

	if ( ! tf_set_ ) {
		loops_operation_->set_design_loop(false);
		loops_operation_->set_include_neighbors(include_neighbor_sc);
		loops_operation_->set_loops(loops);
		loops_operation_->set_cutoff_distance(neighbor_dis_);

		tf_->clear();
		tf_->push_back(cmd_line_operation_ );
		tf_->push_back( restrict_design_operation_ );
		tf_->push_back(loops_operation_);
		tf_->push_back( utility::pointer::make_shared< ExtraRotamers >( 0 /*all*/, 1 /*ex1*/, 1 /*level*/ ) );
		tf_->push_back( utility::pointer::make_shared< ExtraRotamers >( 0 /*all*/, 2 /*ex2*/, 1 /*level*/ ) );
		tf_->push_back(utility::pointer::make_shared< IncludeCurrent >());

	} else {
		tf_->push_back( restrict_design_operation_ );
	}


	BackrubProtocolOP backrub_protocol = utility::pointer::make_shared< BackrubProtocol >();
	backrub_protocol->set_movemap(mm);
	backrub_protocol->set_scorefunction(local_scorefxn);
	backrub_protocol->set_taskfactory(tf_);

	backrub_protocol->apply(pose);
	//protocols::loops::remove_cutpoint_variants(pose);
	//pose.fold_tree(original_ft);

	//local_scorefxn->show(TR, pose);
	TR << "CDR backrub start: "<<start_e << std::endl;
	TR << "CDR backrub end: "<<(*scorefxn_)(pose) << std::endl;
}

void
GeneralAntibodyModeler::repack_antigen_ab_interface(Pose& pose) const {
	vector1< std::string > antigen_chains = ab_info_->get_antigen_chains();
	if ( antigen_chains.size() == 0 ) {
		TR <<" Antigen not present to repack interface" << std::endl;
		return;
	}

	core::kinematics::FoldTree original_ft = pose.fold_tree();
	apply_A_LH_foldtree(pose);

	core::Real start_e = (*scorefxn_)(pose);
	core::pack::task::PackerTaskOP task;
	if ( !tf_set_ ) {
		task = interface_tf_->create_task_and_apply_taskoperations(pose);
	} else {
		task = tf_->create_task_and_apply_taskoperations(pose);
	}

	PackRotamersMover packer = PackRotamersMover(scorefxn_, task);
	packer.apply(pose);
	TR << "repack antigen-antibody interface start: "<<start_e << std::endl;
	TR << "repack antigen-antibody interface end: "<<(*scorefxn_)(pose) << std::endl;
	pose.fold_tree(original_ft);
}

void
GeneralAntibodyModeler::repack_antigen_interface(Pose & pose) const {
	vector1< std::string > antigen_chains = ab_info_->get_antigen_chains();
	if ( antigen_chains.size() == 0 ) {
		TR <<" Antigen not present to repack interface" << std::endl;
		return;
	}

	core::kinematics::FoldTree original_ft = pose.fold_tree();
	apply_A_LH_foldtree(pose);

	core::Real start_e = (*scorefxn_)(pose);
	core::pack::task::TaskFactoryOP tf = interface_tf_->clone();
	core::Size L_chain  = core::pose::get_chain_id_from_chain("L", pose);
	core::Size H_chain = core::pose::get_chain_id_from_chain("H", pose);
	tf->push_back(utility::pointer::make_shared< PreventChainFromRepackingOperation >(L_chain));
	tf->push_back(utility::pointer::make_shared< PreventChainFromRepackingOperation >(H_chain));

	core::pack::task::PackerTaskOP task = tf->create_task_and_apply_taskoperations(pose);

	PackRotamersMover packer = PackRotamersMover(scorefxn_, task);
	packer.apply(pose);
	TR << "repack antigen interface start: "<< start_e << std::endl;
	TR << "repack antigen interface end: "<<(*scorefxn_)(pose) << std::endl;
	pose.fold_tree(original_ft);
}

void
GeneralAntibodyModeler::repack_antibody_interface(Pose & pose) const {
	vector1< std::string > antigen_chains = ab_info_->get_antigen_chains();
	if ( antigen_chains.size() == 0 ) {
		TR <<" Antigen not present to detect antibody interface" << std::endl;
		return;
	}

	TR << "Repacking antibody part of interface" << std::endl;
	core::kinematics::FoldTree original_ft = pose.fold_tree();
	apply_A_LH_foldtree(pose);

	core::Real start_e = (*scorefxn_)(pose) ;
	core::pack::task::TaskFactoryOP tf = interface_tf_->clone();

	for ( core::Size i = 1; i <= antigen_chains.size(); ++i ) {
		core::Size chain = core::pose::get_chain_id_from_chain(antigen_chains[i], pose);
		tf->push_back(utility::pointer::make_shared< PreventChainFromRepackingOperation >(chain));
	}
	core::pack::task::PackerTaskOP task = tf->create_task_and_apply_taskoperations(pose);

	PackRotamersMover packer = PackRotamersMover(scorefxn_, task);
	packer.apply(pose);

	TR << "start: "<<start_e<< std::endl;
	TR << "end: "<<(*scorefxn_)(pose) << std::endl;
	pose.fold_tree(original_ft);
}


void
GeneralAntibodyModeler::dock_low_res(Pose& pose, bool pack_interface ) const {
	core::pose::DockingPartners dock_chains = get_dock_chains_from_ab_dock_chains(ab_info_, ab_dock_chains_);
	TR << "Docking " <<dock_chains << " low res " << std::endl;
	vector1< int > movable_jumps(1, 1);
	core::kinematics::FoldTree original_ft = pose.fold_tree();
	//TR<< pose << std::endl;
	//TR << original_ft << std::endl;

	//vector1<char> chains;
	//for (core::Size i = 1; i <= pose.size(); ++i){
	// chains.push_back(core::pose::get_chain_from_chain_id(pose.chain(i), pose));
	//}

	//std::string chains_s(chains.begin(), chains.end());
	//TR << "Chains: " << chains_s << std::endl;
	protocols::docking::setup_foldtree(pose, dock_chains, movable_jumps);

	core::Real start_e = (*scorefxn_)(pose);
	protocols::docking::DockingLowResOP docker( new protocols::docking::DockingLowRes(docking_scorefxn_low_) );
	protocols::simple_moves::SwitchResidueTypeSetMover cen_switch = protocols::simple_moves::SwitchResidueTypeSetMover("centroid");
	protocols::simple_moves::ReturnSidechainMover return_sc = protocols::simple_moves::ReturnSidechainMover(pose);

	cen_switch.apply(pose);
	docker->apply(pose);
	return_sc.apply(pose);

	pose.fold_tree(original_ft);

	if ( pack_interface ) {
		repack_antigen_ab_interface(pose);
	}
	TR << "dock low-res start: "<<start_e<< std::endl;
	TR << "dock low-res end: "<<(*scorefxn_)(pose) << std::endl;
}

void
GeneralAntibodyModeler::dock_high_res(Pose& pose, core::Size first_cycle, core::Size second_cycle) const {

	core::pose::DockingPartners dock_chains = get_dock_chains_from_ab_dock_chains(ab_info_, ab_dock_chains_);
	TR << "Docking " <<dock_chains << " high res " << std::endl;
	vector1< int > movable_jumps(1, 1);
	core::kinematics::FoldTree original_ft = pose.fold_tree();
	protocols::docking::setup_foldtree(pose, dock_chains, movable_jumps);
	core::Real start_e = (*scorefxn_)(pose);
	protocols::docking::DockMCMProtocolOP docker( new protocols::docking::DockMCMProtocol(1, docking_scorefxn_high_, scorefxn_) );

	if ( tf_set_ ) {
		docker->set_task_factory(tf_);
		docker->set_ignore_default_task(true);
		//tf_->create_task_and_apply_taskoperations(pose)->show(TR);
	}
	docker->set_first_cycle(first_cycle);
	docker->set_second_cycle(second_cycle);

	docker->apply(pose);
	TR << "dock high-res start: " << start_e << std::endl;
	TR << "dock high-res end: "<<(*scorefxn_)(pose) << std::endl;
	pose.fold_tree(original_ft);
}

}//design
}//antibody
}//protocols
