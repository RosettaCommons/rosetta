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
/// @details
///  This contains the functions that create initial positions for docking
///  You can either randomize partner 1 or partner 2, spin partner 2, or
///  perform a simple perturbation.
///  Also contains docking mcm protocol
/// @author Monica Berrondo
/// @author Modified by Sergey Lyskov
/// @author Modified by Sid Chaudhury
/// @author Modified by Jacob Corn

#include <protocols/docking/DockMCMCycle.hh>

// Rosetta Headers
#include <core/kinematics/MoveMap.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/moves/JumpOutMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/docking/SidechainMinMover.hh>

//for resfile reading

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// C++ Headers
#include <string>

#include <basic/Tracer.hh>
using basic::T;

#include <protocols/toolbox/task_operations/InterfaceTaskOperation.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>


using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.docking.DockMCMCycle" );

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
	core::scoring::ScoreFunctionOP scorefxn
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
	core::scoring::ScoreFunctionOP scorefxn,
	core::scoring::ScoreFunctionOP scorefxn_pack
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
	core::scoring::ScoreFunctionOP scorefxn,
	core::scoring::ScoreFunctionOP scorefxn_pack
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
	return protocols::moves::MoverOP( new DockMCMCycle(*this) );
}

void
DockMCMCycle::set_task_factory( core::pack::task::TaskFactoryCOP tf )
{
	tf_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory( *tf ) );
}

void DockMCMCycle::set_default()
{

	using namespace basic::options; //quick hack by rhiju
	using namespace basic::options::OptionKeys::docking; // quick hack by rhiju -- later feed this in through dockingprotocol

	// set translation & rotation move magnitudes
	if ( option[ OptionKeys::docking::dock_mcm_trans_magnitude ].user() ) {
		trans_magnitude_ = option[ OptionKeys::docking::dock_mcm_trans_magnitude ]();
	} else {
		trans_magnitude_ = 0.1;
	}

	if ( option[ OptionKeys::docking::dock_mcm_rot_magnitude ].user() ) {
		rot_magnitude_ = option[ OptionKeys::docking::dock_mcm_rot_magnitude ]();
	} else {
		rot_magnitude_ = 5.0;
	}

	rtmin_ = false;
	scmin_ = false;

	// setup scoring with defaults
	if ( scorefxn_ == NULL ) {
		scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
		scorefxn_pack_ = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	}

	// setup the movemap
	movemap_ = core::kinematics::MoveMapOP( new kinematics::MoveMap() );
	movemap_->set_chi( false );
	movemap_->set_bb( false );
	for ( DockJumps::const_iterator it = movable_jumps_.begin(); it != movable_jumps_.end(); ++it ) {
		movemap_->set_jump( *it, true );
	}

	// perhaps call this dock_minimize_bb_res or something.
	if ( option[ bb_min_res ].user() ) {
		utility::vector1< Size > const & min_res = option[ bb_min_res ]();
		for ( Size n = 1; n <= min_res.size(); n++ ) movemap_->set_bb( min_res[n], true );
	}
	if ( option[ sc_min_res ].user() ) {
		utility::vector1< Size > const & min_res = option[ sc_min_res ]();
		for ( Size n = 1; n <= min_res.size(); n++ ) movemap_->set_chi( min_res[n], true );
	}


	// setup values for minimization
	min_tolerance_ = 0.01;
	min_type_ = std::string( "dfpmin_armijo_nonmonotone" );
	nb_list_ = true;

	// setup the mc object
	mc_ = moves::MonteCarloOP( new moves::MonteCarlo( *scorefxn_, 0.8 ) );
	tf_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory );

	// packing information
	repack_period_ = 8;
}

void DockMCMCycle::set_move_map( core::kinematics::MoveMapOP movemap ) { movemap_ = movemap; }
////////////////////////////////////////////////////////////////////////////////
/// @brief
/// @details
///  decides what to call according to options
void DockMCMCycle::apply( core::pose::Pose & pose )
{
	TR.Debug << "in DockMCMCycle.apply" << std::endl;
	// protocols::moves::PyMolMover pymol;   //JQX: comment out, pymolmover will use the random number, not good for debug

	// only set up on first call to the mover
	if ( ! dock_mcm_cycle_ ) {
		init_mc(pose);
	}

	dock_mcm_cycle_->apply( pose );
	// pymol.apply( pose );   //JQX: comment out
	// pymol.send_energy( pose ); //JQX: comment out
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
	/// It is possible for the MonteCarlo instance to be shared among several movers
	/// in other protocols and we do not want to overwrite the history of that instance,
	/// but we must ensure the MonteCarlo instance is ready for use.
	if ( mc_->last_accepted_pose().empty() ) {
		(*scorefxn_)(pose);
		mc_->reset( pose );
	}
	setup_protocol( pose );
	TR << "Setting up defaults for DockMCMCycle: " << std::endl;
	show(TR);
}

///////////////////////////////////////////////////////////////////////////////////
///
/// @brief main entrance to do monte carlo minimization
/// @details
///   single cycle of mcm minimization.  Then it is rigid-body minimized
///   to a stringent tolerance.
///
/// @remarks
///
/// @references docking_mcm_protocol from docking_minimize.cc
///    pose_docking_monte_carlo_minimize from pose_docking.cc
///
/// @author Monica Berrondo
///
/////////////////////////////////////////////////////////////////////////////////
void DockMCMCycle::setup_protocol( core::pose::Pose & pose ) {
	using namespace moves;
	using namespace basic::options;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace protocols::toolbox::task_operations;

	//JQX: set up rigid body movers
	rigid::RigidBodyPerturbMoverOP rb_mover( new rigid::RigidBodyPerturbMover( pose, *movemap_, rot_magnitude_, trans_magnitude_ , rigid::partner_downstream, true ) );


	//set up sidechain movers for each movable jump
	// tf_->push_back( new RestrictToInterface( movable_jumps_ ) );
	// JQX commented this out, there is one more RestrictToInterface function in the DockTaskFactory.cc file


	protocols::simple_moves::RotamerTrialsMoverOP rottrial( new protocols::simple_moves::RotamerTrialsMover( scorefxn_pack_, tf_ ) );
	SequenceMoverOP rb_pack_min( new SequenceMover );

	rb_pack_min->add_mover( rb_mover );
	rb_pack_min->add_mover( rottrial );


	//JQX: use   (SequenceMover) rb_pack_min   and   (MinMover) min_mover
	//JQX: to define the JumpOutMover
	core::Real minimization_threshold = 15.0;
	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( movemap_, scorefxn_, min_type_, min_tolerance_, nb_list_ ) );
	JumpOutMoverOP rb_mover_min( new JumpOutMover( rb_pack_min, min_mover, scorefxn_, minimization_threshold ) );


	//JQX: (JumpOutMover) rb_mover_min    should be a TrialMover
	TrialMoverOP rb_mover_min_trail( new TrialMover( rb_mover_min, mc_) );


	//JQX: define the SequenceMover repack_step
	//JQX: repack_step is almost the same as rb_mover_min_trail
	//JQX: the only difference is you have one more step: (TrialMover) the pack_interface_and_loops_trial
	//JQX: the TrialMover is actually    ...       (PackRotamersMover) pack_rotamers
	SequenceMoverOP repack_step( new SequenceMover );
	repack_step->add_mover(rb_mover_min_trail);
	protocols::simple_moves::PackRotamersMoverOP pack_rotamers( new protocols::simple_moves::PackRotamersMover( scorefxn_pack_ ) );
	pack_rotamers->task_factory(tf_);
	TrialMoverOP pack_interface_and_move_loops_trial( new TrialMover( pack_rotamers, mc_) );
	repack_step->add_mover(pack_interface_and_move_loops_trial);




	//  these are not being used at all in the extreme code week, JQX incorporated into the sequence
	if ( rtmin_ ) {
		protocols::simple_moves::RotamerTrialsMinMoverOP rtmin( new protocols::simple_moves::RotamerTrialsMinMover( scorefxn_pack_, tf_ ) );
		TrialMoverOP rtmin_trial( new TrialMover( rtmin, mc_ ) );
		repack_step->add_mover(rtmin_trial);
	}
	if ( scmin_ ) {
		core::pack::task::TaskFactoryCOP my_tf( tf_);
		//@TODO JQX: this is so weird, I cannot directly put tf_ to construct the SideChainMinMover
		protocols::docking::SidechainMinMoverOP scmin_mover( new protocols::docking::SidechainMinMover(scorefxn_pack_,  my_tf  ) );
		TrialMoverOP scmin_trial( new TrialMover( scmin_mover, mc_ ) );
		repack_step->add_mover(scmin_trial);
	}

	//JQX: define the cycle mover
	//JQX: 1. rb_mover_min_trail (7 times)
	//JQX: 2. repack_tep (1 time)
	dock_mcm_cycle_ = moves::CycleMoverOP( new CycleMover );
	for ( Size i=1; i<repack_period_; ++i ) dock_mcm_cycle_->add_mover( rb_mover_min_trail );
	dock_mcm_cycle_->add_mover( repack_step );


}

/// @details  Show the complete setup of the docking protocol
void
DockMCMCycle::show( std::ostream & out ) const  {
	out << *this;
}

std::ostream & operator<<(std::ostream& out, const DockMCMCycle & dp )
{
	using namespace ObjexxFCL::format;

	// All output will be 80 characters - 80 is a nice number, don't you think?
	std::string line_marker = "///";
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	// Display the movable jumps that will be used in docking
	out << line_marker << " Dockable Jumps: ";

	core::Size spaces_so_far = 23;
	bool first = true;
	for ( DockJumps::const_iterator it = dp.movable_jumps_.begin() ; it != dp.movable_jumps_.end() ; ++it ) {
		if ( !first ) {
			out << ", ";
			spaces_so_far += 2;
		} else first = false;

		out << I( 1, *it );
		spaces_so_far += 1;
	}
	core::Size remaining_spaces = 80 - spaces_so_far;
	if ( remaining_spaces > 0 ) out << space( 80 - spaces_so_far );
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

DockJumps const & DockMCMCycle::movable_jumps() const
{
	return movable_jumps_;
}

} // namespace docking
} // namespace protocols
