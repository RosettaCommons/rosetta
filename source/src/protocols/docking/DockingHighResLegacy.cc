// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DockingHighRes
/// @brief protocols that are specific to high resolution docking
/// @detailed
///		This contains the functions that create initial positions for docking
///		You can either randomize partner 1 or partner 2, spin partner 2, or
///		perform a simple perturbation.
/// 	Also contains docking mcm protocol
/// @author Monica Berrondo
/// @author Modified by Sergey Lyskov
/// @author Modified by Sid Chaudhury
/// @author Modified by Jacob Corn

#include <protocols/docking/DockingHighResLegacy.hh>
#include <protocols/docking/SidechainMinMover.hh>

// Rosetta Headers
#include <core/kinematics/MoveMap.hh>

// AUTO-REMOVED #include <core/conformation/Interface.hh>

#include <basic/options/option.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictChainToRepackingOperation.hh>
#include <core/conformation/Residue.hh> // for design() flag
#include <core/pack/task/operation/NoRepackDisulfides.hh>
// AUTO-REMOVED #include <core/pack/task/operation/OperateOnCertainResidues.hh>
// AUTO-REMOVED #include <core/pack/task/operation/ResLvlTaskOperations.hh> // PreventRepackingRLT
// AUTO-REMOVED #include <core/pack/task/operation/ResFilters.hh> // ResidueLacksProperty
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>

#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
// AUTO-REMOVED #include <protocols/moves/OutputMovers.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/JumpOutMover.hh>
#include <protocols/moves/ChangeFoldTreeMover.hh>
#include <protocols/moves/RepeatMover.hh>
//for resfile reading
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/LoopMover.fwd.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_Backrub.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <string>

//Utility Headers
// AUTO-REMOVED #include <utility/tag/Tag.hh> // REQUIRED FOR WINDOWS
// AUTO-REMOVED #include <numeric/conversions.hh>

#include <numeric/trig.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>

#include <basic/Tracer.hh>

using basic::T;

// option key includes

#include <basic/options/keys/docking.OptionKeys.gen.hh>

#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>



using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.docking.DockingHighRes");

//     originally from dock_structure.cc Jeff Gray April 2001

using namespace core;

namespace protocols {
namespace docking {

// default constructor
DockingHighResLegacy::DockingHighResLegacy() : DockingHighRes()
{
	// The following three lines are done by the parent class's constructor
	//movable_jumps_.push_back(1); // operate on the first jump
	//scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" ) ;
	//scorefxn_pack() = core::scoring::getScoreFunctionLegacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	moves::Mover::type( "DockingHighResLegacy" );
	init_task_factory_ = NULL;
	design_ = false;
}

// constructor with arguments
// only one movable jump
DockingHighResLegacy::DockingHighResLegacy(
	int rb_jump,
	core::scoring::ScoreFunctionOP scorefxn

) : DockingHighRes( rb_jump, scorefxn )
{
	//movable_jumps_.push_back( rb_jump_in );
	moves::Mover::type( "DockingHighResLegacy" );
	init_task_factory_ = NULL;
	design_ = false;
}

// constructor with arguments
// only one movable jump, scoring and packing defined
DockingHighResLegacy::DockingHighResLegacy(
	int rb_jump,
	core::scoring::ScoreFunctionOP scorefxn,
	core::scoring::ScoreFunctionOP scorefxn_pack
) : DockingHighRes( rb_jump, scorefxn, scorefxn_pack )
{
	//movable_jumps_.push_back( rb_jump_in );
	moves::Mover::type( "DockingHighResLegacy" );
	init_task_factory_ = NULL;
	design_ = false;
}

// constructor with arguments
// only one movable jump, scoring and packing defined
DockingHighResLegacy::DockingHighResLegacy(
	DockJumps const movable_jumps,
	core::scoring::ScoreFunctionOP scorefxn,
	core::scoring::ScoreFunctionOP scorefxn_pack
) : DockingHighRes( movable_jumps, scorefxn, scorefxn_pack )
{
	//movable_jumps_ = movable_jumps;
	moves::Mover::type( "DockingHighResLegacy" );
	init_task_factory_ = NULL;
	design_ = false;
}

//destructor
DockingHighResLegacy::~DockingHighResLegacy() {}

//clone
protocols::moves::MoverOP DockingHighResLegacy::clone() const {
	return new DockingHighResLegacy(*this);
}

void DockingHighResLegacy::set_min_type( std::string min_type_in ) { min_type_ = min_type_in;}

void DockingHighResLegacy::set_repack( bool repack_switch){ repack_switch_ = repack_switch;}

void DockingHighResLegacy::set_trans_magnitude( core::Real trans_magnitude) { trans_magnitude_ = trans_magnitude;}

void DockingHighResLegacy::set_rot_magnitude( core::Real rot_magnitude) { rot_magnitude_ = rot_magnitude;}

void DockingHighResLegacy::set_task_factory( core::pack::task::TaskFactoryOP task )
{
	init_task_factory_ = task;
}

moves::MonteCarloOP DockingHighResLegacy::get_mc() { return mc_; }

void DockingHighResLegacy::design( bool const des ) {
	design_ = des;
}

bool DockingHighResLegacy::design() const { return design_; }

void DockingHighResLegacy::set_default( core::pose::Pose & pose ) {
	using namespace basic::options;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;


	// SETS UP the stuff in pose
	(*scorefxn())( pose );

	core::Real trans_magnitude=option[ OptionKeys::docking::dock_mcm_trans_magnitude ]();
	set_trans_magnitude(trans_magnitude);
	core::Real rot_magnitude=option[ OptionKeys::docking::dock_mcm_rot_magnitude ]();
	set_rot_magnitude(rot_magnitude);

	temperature_ = option[ OptionKeys::docking::temperature]();
	//temperature_ = 0.8;
	repack_switch_ = true;
	repack_period_ = option[ OptionKeys::docking::repack_period]();
	//repack_period_ = 8;

	//sets up MC object
	mc_ = new moves::MonteCarlo( pose, *scorefxn(), temperature_ );

	//sets up default movemap
	bb_ = false;
	chi_ = false;
	movemap_ = new kinematics::MoveMap();
	movemap_->set_chi( chi_ );
	movemap_->set_bb( bb_ );
	for( DockJumps::const_iterator it = movable_jumps().begin(); it != movable_jumps().end(); ++it ) {
		movemap_->set_jump( *it, true );
	}

	//sets up minimization parameters
	min_tolerance_ = 0.01;
	min_type_ = std::string( "dfpmin_armijo_nonmonotone" );
	nb_list_ = true;

	setup_packing( pose );

	//set up docking protocol based on options
	set_protocol( pose );
}


void DockingHighResLegacy::set_move_map(core::kinematics::MoveMapOP movemap_in){ movemap_ = movemap_in; }

void DockingHighResLegacy::set_protocol( pose::Pose & pose )
{
	using namespace basic::options;

	if ( option[ OptionKeys::docking::dock_min ]() ) {
		set_dock_min_protocol();
	} else if ( option[ OptionKeys::docking::dock_ppk ]() ){
		//set_dock_ppk_protocol( pose );
		utility::exit( "Danger Will Robinson! Prepacking is no longer performed within the DockingProtocol.  Use the docking_prepack_protocol executable first, then use the output as the starting structure in docking.", 1 );
	} else {
		set_dock_mcm_protocol( pose );
	}
}

void DockingHighResLegacy::define_loops( pose::Pose const & pose, loops::LoopsOP loop_set, Real & interface_dist ) {
	//runtime_assert( movable_jumps_.size() == 1 ); // CURRENTLY ONLY SUPPORTED WITH SIMPLE DOCKING
	//core::Size const rb_jump = movable_jumps_[1];

    TR << "This method clears the loop set the is passed in.  Before using this, please investigate it thoroughly." << std::endl;

	loop_set->clear();

	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace protocols::toolbox::task_operations;
	//RestrictTaskForDockingOP rtfd = new RestrictTaskForDocking( scorefxn_, rb_jump_, true, interface_dist );
	pack::task::TaskFactory tf;
	//tf.push_back( rtfd );
	//tf.push_back( new RestrictTaskForDocking( scorefxn_, rb_jump_, true, interface_dist ) );
	tf.push_back( new RestrictToInterface( movable_jumps(), interface_dist ) );
	pack::task::PackerTaskOP task = tf.create_task_and_apply_taskoperations( pose );

	// extend one residue beyond borders of repackable regions, don't allow 1-residue loops
	core::Size const nres = pose.total_residue();
	utility::vector1<bool> flexible_region( nres, false );
	for ( Size i=2; i < nres; ++i ) {
		int num_flexible(0);
		if ( pose.pdb_info()->chain(i-1) == pose.pdb_info()->chain(i+1) ) {
			if ( task->pack_residue(i-1) ) ++num_flexible;
			if ( task->pack_residue(i) ) ++num_flexible;
			if ( task->pack_residue(i+1) ) ++num_flexible;
		}
		if ( num_flexible > 1 ) {
			flexible_region.at(i-1) = true;
			flexible_region.at(i) = true;
			flexible_region.at(i+1) = true;
		}
	}

	// if we have a single fixed residue between two loops, make this flexible
	for ( Size i=2; i < nres; ++i ) {
		if ( flexible_region.at(i-1) && flexible_region.at(i+1) ) flexible_region.at(i) = true;
	}

	// jk For now, don't let the first or last two residues of a chain be flexible
	flexible_region.at(1) = false;
	flexible_region.at(2) = false;
	for ( Size i=3; i < (nres-1); ++i ) {
		if ( pose.pdb_info()->chain(i-1) != pose.pdb_info()->chain(i) ) {
			flexible_region.at(i-2) = false;
			flexible_region.at(i-1) = false;
			flexible_region.at(i) = false;
			flexible_region.at(i+1) = false;
		}
	}
	flexible_region.at(nres-1) = false;
	flexible_region.at(nres) = false;

	// disallow one-residue loops
	for ( Size i=2; i < nres; ++i ) {
		if ( ( ! flexible_region.at(i-1)) && ( ! flexible_region.at(i+1)) ) flexible_region.at(i) = false;
	}
	// disallow two-residue loops
	for ( Size i=3; i < nres; ++i ) {
		if ( ( ! flexible_region.at(i-2)) && ( ! flexible_region.at(i+1)) ) {
			flexible_region.at(i-1) = false;
			flexible_region.at(i) = false;
		}
	}

	// setup loops
	core::Size loop_start=0;
	core::Size loop_stop=0;
	for ( Size i=1; i < nres; ++i ) {
		if ( flexible_region.at(i) ) {
			loop_start = i;
			loop_stop = i;
			for ( Size j=i+1; j <= nres; ++j ) {
				// if j is on a different chain than i, break
				if ( pose.pdb_info()->chain(j) != pose.pdb_info()->chain(i) ) {
					break;
				}
				// if j is not flexible, break
				if ( ( ! flexible_region.at(j) ) ) {
					break;
				}
				loop_stop = j;
			}
			loop_set->add_loop( loop_start, loop_stop, 0 );
			i = loop_stop;
		}
	}

	loop_set->choose_cutpoints( pose );

	return;
}


////////////////////////////////////////////////////////////////////////////////
/// @begin docking high resolution apply function
/// @brief
/// @detailed
///		decides what to call according to options
void DockingHighResLegacy::apply( core::pose::Pose & pose )
{
	using namespace scoring;
	using namespace basic::options;

	TR << "in DockingHighRes.apply" << std::endl;

	// jec sanity check to avoid overwriting newly-set minimizers on every apply
	// Need to call set_default to reset the cycle mover index in the case where a decoy fails the filters.
	//if ( !mc_ ) {
		set_default( pose );
	//}

	mc_->reset( pose );

	docking_highres_protocol_mover_->apply( pose );

//	pose.dump_pdb("test.pdb");  // JQX: for testing purpose
//	exit(-1);   // JQX: testing

	if ( !option[ OptionKeys::docking::dock_ppk ]() ) mc_->recover_low( pose );
}

///////////////////////////////////////////////////////////////////////////////////
/// @begin minimize_trial
///
/// @brief main entrance for normal rigid-body minimization
/// @detailed
///		retrieve the structure in the low array and do the normal minimization
///		by calling using a min_mover to optimize the score accourding to the
///		scorefunction that has been set
///
/// @remarks
///
/// @references docking_minimize_trial from docking_minimize.cc
///				pose_docking_minimize_trial from pose_docking.cc
///
/// @authors Monica Berrondo June 14 2007
///
/// @last_modified October 15 2007
/////////////////////////////////////////////////////////////////////////////////
void DockingHighResLegacy::set_dock_min_protocol() {
	using namespace moves;

	TR << "::::::::::::::::::DOCK_MIN:::::::::::::::::::" << std::endl;

	protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( movemap_, scorefxn(), min_type_, min_tolerance_, nb_list_ );
	TrialMoverOP minimize_trial = new TrialMover( min_mover, mc_ );
	docking_highres_protocol_mover_ = new SequenceMover;
	docking_highres_protocol_mover_->add_mover( minimize_trial );
}

///////////////////////////////////////////////////////////////////////////////////
/// @begin dock_mcm_protocol
///
/// @brief main entrance to do monte carlo minimization
/// @detailed
///			a total of 50 cycles of monte-carlo minimization will be
///			carried out if the minimized structure can pass the filter
///			after the first and fifth cycle.  Then it is rigid-body minimized
///			to a stringent tolerance.
///
/// @remarks
///
/// @references docking_mcm_protocol from docking_minimize.cc
///				pose_docking_monte_carlo_minimize from pose_docking.cc
///
/// @authors Sid Chaudhury May 28 2009
///
/// @last_modified April 30 2008
/////////////////////////////////////////////////////////////////////////////////
void DockingHighResLegacy::set_dock_mcm_protocol( core::pose::Pose & pose ) {
	using namespace moves;
	using namespace basic::options;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace protocols::toolbox::task_operations;

	//set up minimizer movers
	protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( movemap_, scorefxn(), min_type_, min_tolerance_, nb_list_ );

	//set up rigid body movers
	rigid::RigidBodyPerturbMoverOP rb_perturb = new rigid::RigidBodyPerturbMover( pose, *movemap_, rot_magnitude_, trans_magnitude_ , rigid::partner_downstream, true );

	//set up sidechain movers for each movable jump

	TaskFactoryOP tf = new TaskFactory( *task_factory() );
	tf->push_back( new RestrictToInterface( movable_jumps() ) );
	set_task_factory( tf );
	//RotamerTrialsMoverOP pack_rottrial = new protocols::simple_moves::RotamerTrialsMover( scorefxn_pack(), task_factory() );

 	//need to explicitly convert to COP from OP in order for SidechainMinMover to work
 	//besides, need to use init_task_factory_, which is otherwise not used from task_factory()
 	core::pack::task::TaskFactoryCOP ctf=init_task_factory_;

	protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial = new protocols::simple_moves::RotamerTrialsMover( scorefxn_pack(), ctf );
	SequenceMoverOP interface_repack_and_move_loops = new moves::SequenceMover;

	std::string const flex_bb_docking_type = option[ OptionKeys::docking::flexible_bb_docking ]();
	if ( flex_bb_docking_type == "fixedbb" ) {
		// Call pack_rotamers, no backbone movement
		protocols::simple_moves::PackRotamersMoverOP pack_interface_repack = new protocols::simple_moves::PackRotamersMover( scorefxn_pack() );
		//pack_interface_repack->task_factory( task_factory() );
		pack_interface_repack->task_factory( ctf );
		interface_repack_and_move_loops->add_mover( pack_interface_repack );
	} else {

		// Call pack_rotamer before and after loop movement
		protocols::simple_moves::PackRotamersMoverOP pack_interface_repack = new protocols::simple_moves::PackRotamersMover( scorefxn_pack() );
		//pack_interface_repack->task_factory( task_factory() );
		pack_interface_repack->task_factory( ctf );
		interface_repack_and_move_loops->add_mover( pack_interface_repack );

		if ( ( flex_bb_docking_type == "ccd" ) || ( flex_bb_docking_type == "kic" ) ||  ( flex_bb_docking_type == "backrub" ) ) {

			core::kinematics::FoldTree docking_fold_tree( pose.fold_tree() );
			core::kinematics::FoldTree loop_fold_tree = docking_fold_tree;

			// jk Create a copy because we can't modify input pose, and we need energies and the docking fold tree
			// Note: we need the docking fold tree so that we identify the correct "interface" when picking loops (ie. jump #1)
			// jk For now, define loops based on the input structure
			// (so that they're always the same, and we don't have to worry about relaxing a priori, akin to prepacking)
			pose::Pose pose_for_loop_defn = *get_input_pose();
			(*scorefxn_pack())(pose_for_loop_defn);
			pose_for_loop_defn.fold_tree( docking_fold_tree );

			protocols::loops::LoopsOP loop_set = new protocols::loops::Loops();
			Real interface_dist = option[ OptionKeys::docking::flexible_bb_docking_interface_dist ];
			define_loops( pose_for_loop_defn, loop_set, interface_dist );

			loops::loop_mover::LoopMoverOP loop_refine;
			if ( flex_bb_docking_type == "ccd" ) {
				// jk CCD loop refinement (fullatom only)
				TR << "Setting up for ccd loop modeling" << std::endl;
				protocols::loops::fold_tree_from_loops( pose, *loop_set, loop_fold_tree );
				// need to pass a clone of the scorefxn because LoopMover requires a non-const scorefxn
				loop_refine = new loops::loop_mover::refine::LoopMover_Refine_CCD( loop_set, scorefxn_pack()->clone() );
			} else if ( flex_bb_docking_type == "kic" ) {
				// jk KIC loop refinement (fullatom only)
				TR << "Setting up for kinematic (kic) loop modeling" << std::endl;
				protocols::loops::fold_tree_from_loops( pose, *loop_set, loop_fold_tree );
				// need to pass a clone of the scorefxn because LoopMover requires a non-const scorefxn
				loop_refine = new loops::loop_mover::refine::LoopMover_Refine_KIC( loop_set, scorefxn_pack()->clone() );
			} else if ( flex_bb_docking_type == "backrub" ) {
				// jk backrub loop refinement (fullatom only)
				TR << "Setting up for backrub loop modeling" << std::endl;
				// jk backrub can't use segments that span a jump residue, so use a simple fold tree here
				// note: this assumes that the termini are not allowed to move (and in define_loops they aren't)
				loop_fold_tree.simple_tree( pose.total_residue() );
				// need to pass a clone of the scorefxn because LoopMover requires a non-const scorefxn
				loop_refine = new loops::loop_mover::refine::LoopMover_Refine_Backrub( loop_set, scorefxn_pack()->clone() );
			}

			moves::ChangeFoldTreeMoverOP get_loop_ft = new moves::ChangeFoldTreeMover( loop_fold_tree );
			moves::ChangeFoldTreeMoverOP get_docking_ft = new moves::ChangeFoldTreeMover( docking_fold_tree );

			interface_repack_and_move_loops->add_mover( get_loop_ft );
			interface_repack_and_move_loops->add_mover( loop_refine );
			interface_repack_and_move_loops->add_mover( get_docking_ft );

		} else {
			TR << "[ ERROR ] Unknown flexible_bb_docking type: " << flex_bb_docking_type << std::endl;
			exit(1);
		}

		interface_repack_and_move_loops->add_mover( pack_interface_repack );
	}


	TrialMoverOP pack_interface_and_move_loops_trial = new TrialMover( interface_repack_and_move_loops, mc_ );

	protocols::simple_moves::RotamerTrialsMinMoverOP rtmin = new protocols::simple_moves::RotamerTrialsMinMover( scorefxn_pack(), ctf );
	TrialMoverOP rtmin_trial = new TrialMover( rtmin, mc_ );

	//InterfaceSidechainMinMoverOP scmin_mover = new InterfaceSidechainMinMover(rb_jump_, scorefxn_pack() );
	SidechainMinMoverOP scmin_mover = new protocols::docking::SidechainMinMover( scorefxn_pack(), ctf );
	TrialMoverOP scmin_trial = new TrialMover( scmin_mover, mc_ );

	// the standard mcm cycle : rb perturbation->rotamer trials->minimization->MC accept
	SequenceMoverOP rb_mover = new SequenceMover;
	rb_mover->add_mover( rb_perturb );
	if ( repack_switch_ ) rb_mover->add_mover( pack_rottrial );

	SequenceMoverOP min_repack_mover = new SequenceMover;
	min_repack_mover->add_mover( min_mover );
	min_repack_mover->add_mover( pack_rottrial );

	core::Real minimization_threshold = option[ OptionKeys::docking::minimization_threshold]();
	//core::Real minimization_threshold = 15.0;

  JumpOutMoverOP rb_mover_min;
	if ( option[ OptionKeys::docking::extra_rottrial] ) {
		rb_mover_min = new JumpOutMover( rb_mover, min_repack_mover, scorefxn(), minimization_threshold );
	} else {
		rb_mover_min = new JumpOutMover( rb_mover, min_mover, scorefxn(), minimization_threshold );
	}

	TrialMoverOP rb_mover_min_trial = new TrialMover( rb_mover_min, mc_ );

	//every step (rb_mover_min_trial): the standard mcm cycle
	//every 8th step (repack_step): standard mcm cycle + repacking
	//try moving loops too, if desired

	SequenceMoverOP repack_step = new SequenceMover;
	repack_step->add_mover(rb_mover_min_trial);

	if ( repack_switch_ ){
		repack_step->add_mover( pack_interface_and_move_loops_trial );
		if ( rt_min() ) repack_step->add_mover(rtmin_trial);
		if ( sc_min() ) repack_step->add_mover(scmin_trial);
		}

	CycleMoverOP rb_mover_min_trial_repack  = new CycleMover;
	for ( Size i=1; i<repack_period_; ++i ) rb_mover_min_trial_repack->add_mover( rb_mover_min_trial );
	rb_mover_min_trial_repack->add_mover( repack_step );

	//set up initial repack mover
	SequenceMoverOP initial_repack = new SequenceMover;
	initial_repack->add_mover(pack_interface_and_move_loops_trial);
	if ( rt_min() ) initial_repack->add_mover(rtmin_trial);
	if ( sc_min() ) initial_repack->add_mover(scmin_trial);

	//set up initial and final min_trial movers for docking
	TrialMoverOP minimize_trial = new TrialMover( min_mover, mc_ );

	//set up mcm cycles and mcm_repack cycles
	RepeatMoverOP mcm_four_cycles = new RepeatMover( rb_mover_min_trial, option[ OptionKeys::docking::dock_mcm_first_cycles]() );
	RepeatMoverOP mcm_fortyfive_cycles = new RepeatMover( rb_mover_min_trial_repack, option[ OptionKeys::docking::dock_mcm_second_cycles]() );

	//set up protocol mover
	TR << "::::::::::::::::::DOCK_MCM:::::::::::::::::::" << std::endl;

	docking_highres_protocol_mover_ = new SequenceMover;
	if (repack_switch_) docking_highres_protocol_mover_->add_mover( initial_repack );
	docking_highres_protocol_mover_->add_mover( minimize_trial );
	docking_highres_protocol_mover_->add_mover( mcm_four_cycles );
	docking_highres_protocol_mover_->add_mover( mcm_fortyfive_cycles );  // JQX: tempoary comment out, doing test on mcm_four_cycles
	docking_highres_protocol_mover_->add_mover( minimize_trial );        // JQX: tempoary comment out, doing test on mcm_four_cycles
}
///////////////////////////////////////////////////////////////////////////////////
/// @begin dock_ppk_protocol
///
/// @brief main entrance to prepacking for docking
/// @detailed
///			does a full repack of all sidechains for each partner after moving them
///		  away by 500A, then it brings them back together
/// @remarks
///
/// @references 	pose_docking_prepack_protocol from pose_docking.cc
///
/// @authors Sid Chaudhury October 3 2008
///
/// @last_modified October 3 2008
/////////////////////////////////////////////////////////////////////////////////
void DockingHighResLegacy::set_dock_ppk_protocol( core::pose::Pose & pose ) {
	using namespace moves;
	using namespace basic::options;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	TR << "::::::::::::::::::DOCK_PPK:::::::::::::::::::" << std::endl;

	//set up translate-by-axis movers
	Real trans_magnitude = 1000;
	utility::vector1< rigid::RigidBodyTransMoverOP > trans_away_vec;
	for( DockJumps::const_iterator it = movable_jumps().begin(); it != movable_jumps().end(); ++it ) {
		core::Size const rb_jump = *it;
		rigid::RigidBodyTransMoverOP translate_away ( new rigid::RigidBodyTransMover( pose, rb_jump ) );
		translate_away->step_size( trans_magnitude );
		trans_away_vec.push_back( translate_away );
	}

	utility::vector1< rigid::RigidBodyTransMoverOP > trans_back_vec;
	for( DockJumps::const_iterator it = movable_jumps().begin(); it != movable_jumps().end(); ++it ) {
		core::Size const rb_jump = *it;
		rigid::RigidBodyTransMoverOP translate_back ( new rigid::RigidBodyTransMover( pose, rb_jump ) );
		translate_back->step_size( trans_magnitude );
		translate_back->trans_axis().negate();
		trans_back_vec.push_back( translate_back );
	}

	PackerTaskOP task = task_factory()->create_task_and_apply_taskoperations( pose ); // does not include restrict to interface

	protocols::simple_moves::PackRotamersMoverOP prepack_full_repack = new protocols::simple_moves::PackRotamersMover( scorefxn_pack(), task );
	protocols::simple_moves::RotamerTrialsMinMoverOP rtmin_mover = new protocols::simple_moves::RotamerTrialsMinMover( scorefxn_pack(), *task );
	SidechainMinMoverOP scmin_mover = new SidechainMinMover(scorefxn_pack(), task);

	// set up protocol
	docking_highres_protocol_mover_ = new SequenceMover;
	if ( sc_min() ) docking_highres_protocol_mover_->add_mover( scmin_mover );
	for( utility::vector1< rigid::RigidBodyTransMoverOP >::iterator it = trans_away_vec.begin(); it != trans_away_vec.end(); ++it ) {
		docking_highres_protocol_mover_->add_mover( *it );
	}
	docking_highres_protocol_mover_->add_mover( prepack_full_repack );
	if ( rt_min() ) docking_highres_protocol_mover_->add_mover( rtmin_mover );
	if ( sc_min() ) docking_highres_protocol_mover_->add_mover( scmin_mover );
	for( utility::vector1< rigid::RigidBodyTransMoverOP >::iterator it = trans_back_vec.begin(); it != trans_back_vec.end(); ++it ) {
		docking_highres_protocol_mover_->add_mover( *it );
	}
}

// noncost pose for load_unboundrot csts
void DockingHighResLegacy::setup_packing( core::pose::Pose & pose ) {
	using namespace basic::options;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	//set upconstructor packer options
	TaskFactoryOP local_tf = new TaskFactory;

	if ( init_task_factory_ ) {
		TR << "Using user-defined TaskFactory." << std::endl;
		local_tf = new TaskFactory( *init_task_factory_ );
	}

	if( design_ ) {
		TR << "Designing during docking" << std::endl;
		DockJumps repack_chains;
		for( int i = 1; i <= int(pose.num_jump()); ++i ) {
			if( find(movable_jumps().begin(), movable_jumps().end(), i ) == movable_jumps().end() ) { // jump is not movable
				core::Size const cutpoint = pose.fold_tree().cutpoint_by_jump( i );
				core::Size const chain = pose.residue( cutpoint ).chain();
				repack_chains.push_back( chain );
			}
		}
		if( repack_chains.size() > 0 ) {
			for( DockJumps::const_iterator it = repack_chains.begin(); it != repack_chains.end(); ++it ) {
				TR << "Not designing chain " << *it << std::endl;
				local_tf->push_back( new protocols::toolbox::task_operations::RestrictChainToRepackingOperation( *it ) ); //
			}
		}
	}
	else { // default case -- restrict everything to repacking.
		local_tf->push_back( new RestrictToRepacking );
	}
//	tf_->push_back( new OperateOnCertainResidues( new PreventRepackingRLT, new ResidueLacksProperty("PROTEIN") ) );
	local_tf->push_back( new InitializeFromCommandline );
	local_tf->push_back( new IncludeCurrent );
	local_tf->push_back( new NoRepackDisulfides );
	if( option[OptionKeys::packing::resfile].user() ) local_tf->push_back( new ReadResfile );

	// DockingNoRepack only works over the first rb_jump in movable_jumps
	// In a 2-body case this separates 1 & 2 based on the only cutpoint
	// In a multibody case, this separates 1 & 2 based on the first cutpoint
	using namespace protocols::toolbox::task_operations;
	if (option[ OptionKeys::docking::norepack1 ].user()) local_tf->push_back( new DockingNoRepack1( movable_jumps()[1] ) );
	if (option[ OptionKeys::docking::norepack2 ].user()) local_tf->push_back( new DockingNoRepack2( movable_jumps()[1] ) );

	// incorporating Ian's UnboundRotamer operation.
	// note that nothing happens if unboundrot option is inactive!
	core::pack::rotamer_set::UnboundRotamersOperationOP unboundrot = new core::pack::rotamer_set::UnboundRotamersOperation();
	unboundrot->initialize_from_command_line();
	operation::AppendRotamerSetOP unboundrot_operation = new operation::AppendRotamerSet( unboundrot );
	local_tf->push_back( unboundrot_operation );
	core::pack::dunbrack::load_unboundrot(pose); // adds scoring bonuses for the "unbound" rotamers, if any

	// note that RestrictToInterfaceOperation is added during set_dock_mcm_protocol

    set_task_factory( local_tf );
    //DockingHighRes::set_task_factory( local_tf );
}

std::string
DockingHighResLegacy::get_name() const {
	return "DockingHighResLegacy";
}


} // namespace docking
} // namespace protocols
