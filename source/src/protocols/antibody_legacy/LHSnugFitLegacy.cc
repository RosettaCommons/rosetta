// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody/LHSnugFitLegacy.cc
/// @brief Build a homology model of an antibody
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)

#include <protocols/antibody_legacy/LHSnugFitLegacy.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/import_pose/import_pose.hh>

#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/OptH.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
using namespace ObjexxFCL::fmt;

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/docking/SidechainMinMover.hh>
#include <protocols/moves/JumpOutMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/antibody/util.hh>






#include <core/chemical/VariantType.hh>
//JQX:: this header file took care of the "CUTPOINT_LOWER" options below



using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.antibody.LHSnugFitLegacy");




using namespace core;
namespace protocols {
namespace antibody {








// default constructor
LHSnugFitLegacy::LHSnugFitLegacy() : Mover() {

}


LHSnugFitLegacy::LHSnugFitLegacy(loops::LoopsOP loops_in ) : Mover() {
	user_defined_ = true;
	init(loops_in, false);
}


LHSnugFitLegacy::LHSnugFitLegacy(antibody::AntibodyInfoOP antibody_in) : Mover() {
	user_defined_ = true;
	init(antibody_in->get_AllCDRs_in_loopsop(),false);
}

LHSnugFitLegacy::LHSnugFitLegacy(antibody::AntibodyInfoOP antibody_in, bool camelid) : Mover() {
	user_defined_ = true;
	init(antibody_in->get_AllCDRs_in_loopsop(), camelid);
}


// default destructor
LHSnugFitLegacy::~LHSnugFitLegacy() {}

//clone
protocols::moves::MoverOP LHSnugFitLegacy::clone() const {
	return( new LHSnugFitLegacy() );
}





void LHSnugFitLegacy::init(loops::LoopsOP loops_in, bool camelid ) {
	is_camelid_ = camelid;
	all_loops_ = loops_in;
}


void LHSnugFitLegacy::set_default() {
	min_type_="dfpmin_armijo_nonmonotone";
	rot_mag_ = 5.0 ;
	trans_mag_ = 0.1 ;
	temperature_ = 0.8;
}


std::string LHSnugFitLegacy::get_name() const {
	return "LHSnugFitLegacy";
}










void LHSnugFitLegacy::apply( pose::Pose & pose ) {


	using namespace moves;
	bool nb_list = true;
	Size nres = pose.total_residue();

	// rb minimization
	Real min_threshold ( 15.0 ); /* score unit */

	// score functions
	using namespace core::scoring;
	core::scoring::ScoreFunctionOP dock_scorefxn;
	dock_scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
	dock_scorefxn->set_weight( core::scoring::chainbreak, 1.0 );
	dock_scorefxn->set_weight( core::scoring::overlap_chainbreak, 10./3. );

	// score functions
	core::scoring::ScoreFunctionOP pack_scorefxn;
	pack_scorefxn = core::scoring::getScoreFunctionLegacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );

	// remove cutpoints variants for all cdrs
	// "true" forces removal of variants even from non-cutpoints
	loops::remove_cutpoint_variants( pose, true );

	using namespace core::chemical;
	for ( loops::Loops::const_iterator it = all_loops_->begin(),
	        it_end = all_loops_->end();	it != it_end; ++it ) {
		core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, it->cut() );
		core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER,it->cut()+1);
	}

	//setting MoveMap
	kinematics::MoveMapOP cdr_dock_map;
	cdr_dock_map = new kinematics::MoveMap();

	*cdr_dock_map=ab_info_->get_MoveMap_for_LoopsandDock(pose, *ab_info_->get_AllCDRs_in_loopsop(), false, true, 10.0);

	//set up minimizer movers
	simple_moves::MinMoverOP min_mover = new simple_moves::MinMover( cdr_dock_map, dock_scorefxn, min_type_, min_threshold, nb_list );

	//set up rigid body movers
	rigid::RigidBodyPerturbMoverOP rb_perturb = new rigid::RigidBodyPerturbMover( pose,
	        *cdr_dock_map, rot_mag_, trans_mag_, rigid::partner_downstream, true );


	//set up sidechain movers for rigid body jump and loop & neighbors
	utility::vector1_size rb_jump;
	rb_jump.push_back( 1 );
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	// selecting movable c-terminal residues
	utility::vector1< bool> sc_is_flexible( nres, false );
	select_loop_residues( pose, *all_loops_, true/*include_neighbors*/, sc_is_flexible);

	ObjexxFCL::FArray1D_bool loop_residues( nres, false );
	for( Size i = 1; i <= nres; i++ )
		loop_residues( i ) = sc_is_flexible[ i ]; // check mapping
	using namespace protocols::toolbox::task_operations;
	tf_->push_back( new RestrictToInterface( rb_jump, loop_residues ) );



	simple_moves::RotamerTrialsMoverOP pack_rottrial = new simple_moves::RotamerTrialsMover( pack_scorefxn, tf_ );

	simple_moves::PackRotamersMoverOP pack_interface_repack = new simple_moves::PackRotamersMover( pack_scorefxn );
	pack_interface_repack->task_factory(tf_);


	MonteCarloOP mc = new MonteCarlo( pose, *dock_scorefxn, temperature_ );

	TrialMoverOP pack_interface_trial = new TrialMover(pack_interface_repack, mc );

	protocols::docking::SidechainMinMoverOP scmin_mover = new protocols::docking::SidechainMinMover( core::scoring::ScoreFunctionOP( pack_scorefxn ), core::pack::task::TaskFactoryCOP( tf_ ) );
	TrialMoverOP scmin_trial = new TrialMover( scmin_mover, mc );

	SequenceMoverOP rb_mover = new SequenceMover;
	rb_mover->add_mover( rb_perturb );
	rb_mover->add_mover( pack_rottrial );

	JumpOutMoverOP rb_mover_min = new JumpOutMover( rb_mover, min_mover, dock_scorefxn, min_threshold);
	TrialMoverOP rb_mover_min_trial = new TrialMover( rb_mover_min, mc  );

	SequenceMoverOP repack_step = new SequenceMover;
	repack_step->add_mover( rb_mover_min_trial );
	repack_step->add_mover( pack_interface_trial );
	repack_step->add_mover( scmin_trial );

	CycleMoverOP rb_mover_min_trial_repack  = new CycleMover;
	for ( Size i=1; i < 8; ++i )
		rb_mover_min_trial_repack->add_mover( rb_mover_min_trial );
	rb_mover_min_trial_repack->add_mover( repack_step );

	//set up initial repack mover
	SequenceMoverOP initial_repack = new SequenceMover;
	initial_repack->add_mover( pack_interface_trial );
	initial_repack->add_mover( scmin_trial );

	//set up initial and final min_trial movers for docking
	TrialMoverOP minimize_trial = new TrialMover( min_mover, mc );

	//set up mcm cycles and mcm_repack cycles
	RepeatMoverOP mcm_four_cycles = new RepeatMover( rb_mover_min_trial, 4 );

	Size cycles = 3;
	if ( benchmark_ ) cycles = 1;
	RepeatMoverOP mcm_final_cycles = new RepeatMover( rb_mover_min_trial_repack, cycles );

	SequenceMoverOP snugfit_mcm = new SequenceMover;
	snugfit_mcm->add_mover( initial_repack );
	snugfit_mcm->add_mover( minimize_trial );
	snugfit_mcm->add_mover( mcm_four_cycles );
	snugfit_mcm->add_mover( mcm_final_cycles );
	snugfit_mcm->add_mover( minimize_trial );

	snugfit_mcm->apply ( pose );

	return;
}













} // namespace antibody
} // namespace protocols





