// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DockMCMCycle
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

#include <protocols/docking/DockMCMCycle.hh>
// AUTO-REMOVED #include <protocols/docking/SidechainMinMover.hh>
// AUTO-REMOVED #include <protocols/docking/DockTaskFactory.hh>

// Rosetta Headers
#include <core/kinematics/MoveMap.hh>

// AUTO-REMOVED #include <basic/options/option.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <protocols/toolbox/task_operations/RestrictToInterface.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperations.hh>
// AUTO-REMOVED #include <protocols/toolbox/task_operations/RestrictChainToRepackingOperation.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh> // for design() flag
// AUTO-REMOVED #include <core/pack/task/operation/NoRepackDisulfides.hh>
// AUTO-REMOVED #include <core/pack/task/operation/OperateOnCertainResidues.hh>
// AUTO-REMOVED #include <core/pack/task/operation/ResLvlTaskOperations.hh> // PreventRepackingRLT
// AUTO-REMOVED #include <core/pack/task/operation/ResFilters.hh> // ResidueLacksProperty
// AUTO-REMOVED #include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/RotamerConstraint.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/moves/JumpOutMover.hh>
#include <protocols/moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PackRotamersMover.hh>
// AUTO-REMOVED #include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RotamerTrialsMover.hh>
// AUTO-REMOVED #include <protocols/moves/RotamerTrialsMinMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/moves/TrialMover.hh>

//for resfile reading

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// C++ Headers
#include <string>

//Utility Headers
// AUTO-REMOVED #include <utility/tag/Tag.hh> // REQUIRED FOR WINDOWS

#include <basic/Tracer.hh>
using basic::T;

// option key includes
// AUTO-REMOVED #include <basic/options/keys/docking.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <protocols/toolbox/task_operations/InterfaceTaskOperation.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <basic/options/keys/OptionKeys.hh>


using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.docking.DockMCMCycle");

//     originally from dock_structure.cc Jeff Gray April 2001

using namespace core;

namespace protocols {
namespace docking {

// default constructor
DockMCMCycle::DockMCMCycle() : Mover()
{
	movable_jumps_.push_back(1); // operate on the first jump
	moves::Mover::type( "DockMCMCycle" );
	set_default();
}

// constructor with arguments
// only one movable jump
DockMCMCycle::DockMCMCycle(
	core::Size const rb_jump,
	core::scoring::ScoreFunctionCOP scorefxn
) : Mover(), scorefxn_(scorefxn), scorefxn_pack_(scorefxn)
{
	movable_jumps_.push_back( rb_jump );
	moves::Mover::type( "DockMCMCycle" );
	set_default();
}

// constructor with arguments
// only one movable jump, scoring and packing defined
DockMCMCycle::DockMCMCycle(
	core::Size const rb_jump,
	core::scoring::ScoreFunctionCOP scorefxn,
	core::scoring::ScoreFunctionCOP scorefxn_pack
) : Mover(), scorefxn_(scorefxn), scorefxn_pack_(scorefxn_pack)
{
	movable_jumps_.push_back( rb_jump );
	moves::Mover::type( "DockMCMCycle" );
	set_default();
}

// constructor with arguments
// only one movable jump, scoring and packing defined
DockMCMCycle::DockMCMCycle(
	DockJumps const movable_jumps,
	core::scoring::ScoreFunctionCOP scorefxn,
	core::scoring::ScoreFunctionCOP scorefxn_pack
) : Mover(), scorefxn_(scorefxn), scorefxn_pack_(scorefxn_pack)
{
	movable_jumps_ = movable_jumps;
	moves::Mover::type( "DockMCMCycle" );
	set_default();
}

//destructor
DockMCMCycle::~DockMCMCycle() {}

//clone
protocols::moves::MoverOP DockMCMCycle::clone() const {
	return new DockMCMCycle(*this);
}

void
DockMCMCycle::set_task_factory( core::pack::task::TaskFactoryCOP tf )
{
	tf_ = new core::pack::task::TaskFactory( *tf );
}

void DockMCMCycle::set_default()
{

	trans_magnitude_ = 0.1;
	rot_magnitude_ = 5.0;

	// setup scoring with defaults
	if ( scorefxn_() == NULL ) {
		scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
		scorefxn_pack_ = core::scoring::ScoreFunctionFactory::create_score_function( "standard" );
	}

	// setup the movemap
	movemap_ = new kinematics::MoveMap();
	movemap_->set_chi( false );
	movemap_->set_bb( false );
	for ( DockJumps::const_iterator it = movable_jumps_.begin(); it != movable_jumps_.end(); ++it ) {
		movemap_->set_jump( *it, true );
	}

	// setup values for minimization
	min_tolerance_ = 0.01;
	min_type_ = std::string( "dfpmin_armijo_nonmonotone" );
	nb_list_ = true;

	// setup the mc object
	mc_ = new moves::MonteCarlo( *scorefxn_, 0.8 );
	tf_ = new core::pack::task::TaskFactory;

	// packing information
	repack_period_ = 8;
}

void DockMCMCycle::set_move_map( core::kinematics::MoveMapOP movemap ) { movemap_ = movemap; }
////////////////////////////////////////////////////////////////////////////////
/// @begin docking high resolution apply function
/// @brief
/// @detailed
///		decides what to call according to options
void DockMCMCycle::apply( core::pose::Pose & pose )
{
	TR.Debug << "in DockMCMCycle.apply" << std::endl;
//	protocols::moves::PyMolMover pymol;   //JQX: comment out, pymolmover will use the random number, not good for debug

	// only set up on first call to the mover
	if ( mc_->last_accepted_pose().empty() ) {
		init_mc(pose);
	}

	dock_mcm_cycle_->apply( pose );
//	pymol.apply( pose );   //JQX: comment out
//	pymol.send_energy( pose ); //JQX: comment out
}

std::string
DockMCMCycle::get_name() const {
	return "DockMCMCycle";
}




//JQX: reset the pack_cycle_ index
void DockMCMCycle::reset_cycle_index()
{
	dock_mcm_cycle_->reset_cycle_index();
}

void DockMCMCycle::init_mc(core::pose::Pose & pose)
{
	(*scorefxn_)(pose);
	mc_->reset( pose );
	setup_protocol( pose );
	TR << "Setting up defaults for DockMCMCycle: " << std::endl;
	show(TR);
}







///////////////////////////////////////////////////////////////////////////////////
/// @begin dock_mcm_protocol
///
/// @brief main entrance to do monte carlo minimization
/// @detailed
///			single cycle of mcm minimization.  Then it is rigid-body minimized
///			to a stringent tolerance.
///
/// @remarks
///
/// @references docking_mcm_protocol from docking_minimize.cc
///				pose_docking_monte_carlo_minimize from pose_docking.cc
///
/// @authors Monica Berrondo
///
/// @last_modified August 19 2010
/////////////////////////////////////////////////////////////////////////////////
//void DockMCMCycle::setup_protocol( core::pose::Pose & pose ) {
//	using namespace moves;
//	using namespace basic::options;
//	using namespace core::pack::task;
//	using namespace core::pack::task::operation;
//	using namespace protocols::toolbox::task_operations;
//
//	// @TODO these are not being used at all, need to be incorporated into the sequence
////	RotamerTrialsMinMoverOP rtmin = new RotamerTrialsMinMover( scorefxn_pack_, tf_ );
////	TrialMoverOP rtmin_trial = new TrialMover( rtmin, mc_ );
////
////	SidechainMinMoverOP scmin_mover = new SidechainMinMover(scorefxn_pack_, tf_ );
////	TrialMoverOP scmin_trial = new TrialMover( scmin_mover, mc_ );
//
//	//set up rigid body movers
//	rigid::RigidBodyPerturbMoverOP rb_perturb = new rigid::RigidBodyPerturbMover( pose, *movemap_, rot_magnitude_, trans_magnitude_ , rigid::partner_downstream, true );
//
//	//set up sidechain movers for each movable jump
//	tf_->push_back( new RestrictToInterface( movable_jumps_ ) );
//
//	RotamerTrialsMoverOP rottrial = new RotamerTrialsMover( scorefxn_pack_, tf_ );
//
//	// old loop here
//	PackRotamersMoverOP pack_rotamers = new PackRotamersMover( scorefxn_pack_ );
//	pack_rotamers->task_factory(tf_);
//	TrialMoverOP pack_trial = new TrialMover( pack_rotamers, mc_ );
//
//	SequenceMoverOP rb_mover = new SequenceMover;
//	rb_mover->add_mover( rb_perturb );
//	rb_mover->add_mover( rottrial );
//
//	//set up minimizer movers
//	moves::MinMoverOP min_mover = new moves::MinMover( movemap_, scorefxn_, min_type_, min_tolerance_, nb_list_ );
//	core::Real minimization_threshold = 15.0;
//	JumpOutMoverOP rb_min = new JumpOutMover( rb_mover, min_mover, scorefxn_, minimization_threshold );
//	TrialMoverOP rb_min_trial = new TrialMover( rb_min, mc_ );
//
//	SequenceMoverOP rb_min_pack = new SequenceMover;
//	rb_min_pack->add_mover( rb_min_trial );
//	rb_min_pack->add_mover( pack_trial );
//
//	CycleMoverOP pack_cycle = new CycleMover;
//	for (Size i=1; i<repack_period_; ++i) pack_cycle->add_mover( rb_min_trial );
//	pack_cycle->add_mover( rb_min_pack );
//
//	dock_mcm_mover_ = new TrialMover( pack_cycle, mc_ );
//}







/*    JQX:     completely comment this out and rewrite it !!!!!


void DockMCMCycle::setup_protocol( core::pose::Pose & pose ) {
	using namespace moves;
	using namespace basic::options;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
  using namespace protocols::toolbox::task_operations;

	// @TODO these are not being used at all, need to be incorporated into the sequence
//	RotamerTrialsMinMoverOP rtmin = new RotamerTrialsMinMover( scorefxn_pack_, tf_ );
//	TrialMoverOP rtmin_trial = new TrialMover( rtmin, mc_ );
//
//	SidechainMinMoverOP scmin_mover = new SidechainMinMover(scorefxn_pack_, tf_ );
//	TrialMoverOP scmin_trial = new TrialMover( scmin_mover, mc_ );

	//set up rigid body movers
	rigid::RigidBodyPerturbMoverOP rb_mover = new rigid::RigidBodyPerturbMover( pose, *movemap_, rot_magnitude_, trans_magnitude_ , rigid::partner_downstream, true );

	//set up sidechain movers for each movable jump
//	tf_->push_back( new RestrictToInterface( movable_jumps_ ) );   //JQX: temporarly commented out this, because I commented out the Legacy code for this for testing

	RotamerTrialsMoverOP rottrial = new RotamerTrialsMover( scorefxn_pack_, tf_ );

	// old loop here
	PackRotamersMoverOP pack_rotamers = new PackRotamersMover( scorefxn_pack_ );
	pack_rotamers->task_factory(tf_);

	CycleMoverOP pack_cycle = new CycleMover;
	for (Size i=1; i<repack_period_; ++i) pack_cycle->add_mover( rottrial );
	pack_cycle->add_mover( pack_rotamers );

	//set up minimizer movers
	moves::MinMoverOP min_mover = new moves::MinMover( movemap_, scorefxn_, min_type_, min_tolerance_, nb_list_ );

	// the standard mcm cycle : rb perturbation->rotamer trials->minimization->MC accept
	SequenceMoverOP rb_pack_min = new SequenceMover;
	rb_pack_min->add_mover( rb_mover );
	rb_pack_min->add_mover( pack_cycle );
//	rb_pack_min->add_mover( min_mover ); // JQX comment this out, it was pulled out from rb_pack_min, so that I can define JumpOutMover, see below

//	dock_mcm_mover_ = new TrialMover( rb_pack_min, mc_ );  //JQX commented this out, pull out the "min_mover" out of "rb_pack_min", redefine a JumpOutMOver to match Legacy
	core::Real minimization_threshold = 15.0;              // JQX add this, in order to define a new   JumpOutMover
	JumpOutMoverOP rb_mover_min = new JumpOutMover( rb_pack_min, min_mover, scorefxn_, minimization_threshold ); // JQX define the new JumpOutMover called rb_mover_min, match the Legacy code
	dock_mcm_mover_ = new TrialMover( rb_mover_min,mc_ );  //JQX redfine the dock_mcm_mover_
}

*/
//JQX: completely comment above things









// JQX rewrote the setup_protocol below

void DockMCMCycle::setup_protocol( core::pose::Pose & pose ) {
	using namespace moves;
	using namespace basic::options;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace protocols::toolbox::task_operations;



	// @TODO these are not being used at all, need to be incorporated into the sequence
	//	RotamerTrialsMinMoverOP rtmin = new RotamerTrialsMinMover( scorefxn_pack_, tf_ );
	//	TrialMoverOP rtmin_trial = new TrialMover( rtmin, mc_ );
	//
	//	SidechainMinMoverOP scmin_mover = new SidechainMinMover(scorefxn_pack_, tf_ );
	//	TrialMoverOP scmin_trial = new TrialMover( scmin_mover, mc_ );




	//JQX: set up rigid body movers
	rigid::RigidBodyPerturbMoverOP rb_mover = new rigid::RigidBodyPerturbMover( pose, *movemap_, rot_magnitude_, trans_magnitude_ , rigid::partner_downstream, true );


	//set up sidechain movers for each movable jump
	//	tf_->push_back( new RestrictToInterface( movable_jumps_ ) );
			// JQX commented this out, there is one more RestrictToInterface function in the DockTaskFactory.cc file


	RotamerTrialsMoverOP rottrial = new RotamerTrialsMover( scorefxn_pack_, tf_ );
	SequenceMoverOP rb_pack_min = new SequenceMover;

	rb_pack_min->add_mover( rb_mover );
	rb_pack_min->add_mover( rottrial );


	//JQX: use   (SequenceMover) rb_pack_min   and   (MinMover) min_mover
	//JQX: to define the JumpOutMover
	core::Real minimization_threshold = 15.0;
	moves::MinMoverOP min_mover = new moves::MinMover( movemap_, scorefxn_, min_type_, min_tolerance_, nb_list_ );
	JumpOutMoverOP rb_mover_min = new JumpOutMover( rb_pack_min, min_mover, scorefxn_, minimization_threshold );


	//JQX: (JumpOutMover) rb_mover_min    should be a TrialMover
	TrialMoverOP rb_mover_min_trail = new TrialMover( rb_mover_min, mc_);


	//JQX: define the SequenceMover repack_step
	//JQX: repack_step is almost the same as rb_mover_min_trail
	//JQX: the only difference is you have one more step: (TrialMover) the pack_interface_and_loops_trial
	//JQX: the TrialMover is actually    ...       (PackRotamersMover) pack_rotamers
	SequenceMoverOP repack_step = new SequenceMover;
	repack_step->add_mover(rb_mover_min_trail);
	PackRotamersMoverOP pack_rotamers = new PackRotamersMover( scorefxn_pack_ ); pack_rotamers->task_factory(tf_);
	TrialMoverOP pack_interface_and_move_loops_trial = new TrialMover( pack_rotamers, mc_);
	repack_step->add_mover(pack_interface_and_move_loops_trial);


	//JQX: define the cycle mover
	//JQX: 1. rb_mover_min_trail (7 times)
	//JQX: 2. repack_tep (1 time)
	dock_mcm_cycle_ = new CycleMover;
	for (Size i=1; i<repack_period_; ++i) dock_mcm_cycle_->add_mover( rb_mover_min_trail );
	dock_mcm_cycle_->add_mover( repack_step );


}














/// @details  Show the complete setup of the docking protocol
void
DockMCMCycle::show( std::ostream & out ) const  {
	out << *this;
}

std::ostream & operator<<(std::ostream& out, const DockMCMCycle & dp )
{
	using namespace ObjexxFCL::fmt;

	// All output will be 80 characters - 80 is a nice number, don't you think?
	std::string line_marker = "///";
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	// Display the movable jumps that will be used in docking
	out << line_marker << " Dockable Jumps: ";

	core::Size spaces_so_far = 23;
	bool first = true;
	for ( DockJumps::const_iterator it = dp.movable_jumps_.begin() ; it != dp.movable_jumps_.end() ; ++it ){
		if (!first) {
			out << ", ";
			spaces_so_far += 2;
		}
		else first = false;

		out << I( 1, *it );
		spaces_so_far += 1;
	}
	core::Size remaining_spaces = 80 - spaces_so_far;
	if ( remaining_spaces > 0 )	out << space( 80 - spaces_so_far );
	out << line_marker << std::endl;

	// Display the translation, and rotation of the cycle
	out << line_marker << " Translation: " << dp.trans_magnitude_<< space( 57 ) << line_marker << std::endl;
	out << line_marker << " Rotation: " << dp.rot_magnitude_ << space( 62 ) << line_marker << std::endl;
	out << line_marker << " Scorefunction: " << space( 58 ) << line_marker << std::endl;
	dp.scorefxn_->show(out);
	out <<std::endl;
	out << line_marker << " Packing scorefunction: " << space( 50 ) << line_marker << std::endl;
	dp.scorefxn_pack_->show(out);
	out << std::endl;

	// Close the box I have drawn
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	return out;
}


} // namespace docking
} // namespace protocols
