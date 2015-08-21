// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DockingLowRes
/// @brief protocols that are specific to docking low resolution
/// @details
/// @author Monica Berrondo
/// @author Modified by Sergey Lyskov
/// @author Modified by Jacob Corn
/// @author Modified by Daisuke Kuroda

#include <protocols/docking/DockingLowRes.hh>

// Rosetta Headers
#include <core/kinematics/MoveMap.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>


#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/scoring/InterfaceInfo.hh>

// Utility Headers
#include <utility/tools/make_vector1.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// C++ Headers
#include <string>

//Utility Headers

#include <numeric/trig.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>

#include <basic/Tracer.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.docking.DockingLowRes" );

//     originally from dock_structure.cc Jeff Gray April 2001
//

using namespace core;

namespace protocols {
namespace docking {

// default constructor
DockingLowRes::DockingLowRes() :
	Mover()
{
	init(utility::tools::make_vector1<core::SSize>(1), NULL);
}


// constructor with arguments
DockingLowRes::DockingLowRes(
	core::scoring::ScoreFunctionCOP scorefxn,
	core::Size const rb_jump
) : Mover(), scorefxn_(scorefxn)
{
	init(utility::tools::make_vector1<core::SSize>(rb_jump), scorefxn);
}

DockingLowRes::DockingLowRes(
	core::scoring::ScoreFunctionCOP scorefxn,
	DockJumps const movable_jumps
) : Mover(), scorefxn_(scorefxn)
{
	init(movable_jumps, scorefxn);
}

void DockingLowRes::init(
	DockJumps const movable_jumps,
	core::scoring::ScoreFunctionCOP scorefxn
)
{
	moves::Mover::type( "DockingLowRes" );
	movable_jumps_ = movable_jumps;

	// setup all the booleans with default values
	// they will get overwritten by the options and/or passed values
	set_default();

	if ( scorefxn == NULL ) {
		scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "interchain_cen" );
	} else scorefxn_ = scorefxn;

	// set up objects based on the boolean values defined above
	sync_objects_with_flags();
}


//destructor
DockingLowRes::~DockingLowRes() {}

protocols::moves::MoverOP
DockingLowRes::clone() const {
	return protocols::moves::MoverOP( new DockingLowRes(*this) );
}

void DockingLowRes::set_trans_magnitude( core::Real trans_magnitude) { trans_magnitude_ = trans_magnitude;}

void DockingLowRes::set_rot_magnitude( core::Real rot_magnitude) { rot_magnitude_ = rot_magnitude;}


void DockingLowRes::set_default() {

	core::Real trans_magnitude=basic::options::option[ basic::options::OptionKeys::docking::docklowres_trans_magnitude ]();
	set_trans_magnitude(trans_magnitude);
	core::Real rot_magnitude=basic::options::option[ basic::options::OptionKeys::docking::docklowres_rot_magnitude ]();
	set_rot_magnitude(rot_magnitude);
	//TR << " Lowres docking trans_magnitude_ " << trans_magnitude_ << std::endl;
	//TR << " Lowres docking rot_magnitude_" << rot_magnitude_ << std::endl;

	//trans_magnitude_ = 0.7;
	//rot_magnitude_ = 5.0;
	chi_ = false;
	bb_ = false;

	temperature_ = 0.8;
	nb_list_ = true; /// not sure if this should be true or not
	accept_rate_ = 0.0;

	inner_cycles_ = 50;
	outer_cycles_ = 10;
}

void DockingLowRes::sync_objects_with_flags()
{
	rb_mover_ = NULL;
	docking_lowres_protocol_ = NULL;

	// the movable dof's -- jumps only in this case
	movemap_ = core::kinematics::MoveMapOP( new kinematics::MoveMap() );
	movemap_->set_chi( chi_ ); // is this right?
	movemap_->set_bb( bb_ ); // is this right?
	for ( DockJumps::const_iterator it = movable_jumps_.begin(); it != movable_jumps_.end(); ++it ) {
		movemap_->set_jump( *it, true );
	}

	// setup the mc object
	mc_ = protocols::moves::MonteCarloOP( new moves::MonteCarlo( *scorefxn_, temperature_ ) );

	flags_and_objects_are_in_sync_ = true;
	first_apply_with_current_setup_ = true;
}

void DockingLowRes::set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn )
{
	scorefxn_ = scorefxn;
	// mc object score function is out of sync and needs to be recreated
	flags_and_objects_are_in_sync_ = false;
}

void DockingLowRes::finalize_setup( core::pose::Pose & pose){
	using namespace moves;

	rb_mover_ = protocols::rigid::RigidBodyPerturbNoCenterMoverOP( new rigid::RigidBodyPerturbNoCenterMover( pose, *movemap_, rot_magnitude_, trans_magnitude_, protocols::rigid::n2c ) );

	docking_lowres_protocol_ = protocols::moves::SequenceMoverOP( new SequenceMover );
	docking_lowres_protocol_->add_mover( rb_mover_ );

	//JAB - moved from setup_foldtree.  It should be here.

	//set up InterfaceInfo object in pose to specify which interface(s)
	//to calculate docking centroid mode scoring components from
	using namespace core::scoring;
	using namespace protocols::scoring;
	//using core::pose::datacache::CacheableDataType::INTERFACE_INFO;

	InterfaceInfoOP docking_interface( new InterfaceInfo( movable_jumps_ ) );
	pose.data().set(
		core::pose::datacache::CacheableDataType::INTERFACE_INFO,
		docking_interface);
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief Perform several cycles of rigid-body Monte Carlo moves
///       and adapt the step size.
/// @details
///
/// @remarks
///       currently used only in the low-resolution step (centroid mode)
///
/// @references pose_docking_centroid_rigid_body_adaptive from pose_docking.cc and
///    rigid_body_MC_cycle_adaptive from dock_structure.cc
///
/// @author Monica Berrondo October 22 2007
///
/////////////////////////////////////////////////////////////////////////////////
void DockingLowRes::apply( core::pose::Pose & pose )
{
	using namespace scoring;

	//TR << "in " << name << "apply" << std::endl;
	TR << "in DockingLowRes.apply" << std::endl;

	if ( !flags_and_objects_are_in_sync_ ) {
		sync_objects_with_flags();
	}

	if ( first_apply_with_current_setup_ ) {
		finalize_setup( pose );
		first_apply_with_current_setup_ = false;
	}

	/// since the mc object is created on construction, the pose must be passed to it
	/// every time the mover is applied to ensure that the mc object is in sync with the pose
	(*scorefxn_)( pose );
	mc_->reset( pose );

	show( TR );

	TR << "::::::::::::::::::Centroid Rigid Body Adaptive:::::::::::::::::::" << std::endl;

	for ( core::Size i=1; i<=outer_cycles_; ++i ) {
		rigid_body_trial( pose );
		if ( accept_rate_ < 0.5 ) {
			trans_magnitude_ *= 0.9;
			rot_magnitude_ *= 0.9;
		} else {
			trans_magnitude_ *= 1.1;
			rot_magnitude_ *= 1.1;
		}
		//  pose.energies().show( std::cout );
	}
	mc_->recover_low( pose );
	TR.flush();
	// pose.energies().show( std::cout );
}

std::string
DockingLowRes::get_name() const {
	return type();
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief Perform a cycle of rigid-body Monte Carlo moves
///
/// @details  Performs a number (nattempts) of MC rigid-body moves
///       (of size trans_magnitude, rot_magnitude). The number of successful
///       attempts is stored in accept_rate_ and used in adaptive trials.
///
/// @remarks the success_rate defines
///       whether the translation/rotation size is increased or decreased for
///       the next cycle.
///       currently used only in the low-resolution step (centroid mode)
///
/// @references pose_docking_rigid_body_trial from pose_docking.cc and
///    rigid_body_MC_cycle from dock_structure.cc
///
/// @author Monica Berrondo October 22 2007
///
/////////////////////////////////////////////////////////////////////////////////
void DockingLowRes::rigid_body_trial( core::pose::Pose & pose )
{
	using namespace moves;

	// PDBDumpMoverOP dump = new PDBDumpMover("lowres_cycle_");
	// dump->apply( pose );
	// MCShowMoverOP mc_show = new MCShowMover( mc_ );
	// mc_show->apply( pose );

	rb_mover_->rot_magnitude( rot_magnitude_ );
	rb_mover_->trans_magnitude( trans_magnitude_ );

	TrialMoverOP rb_trial( new TrialMover( docking_lowres_protocol_, mc_ ) );
	// rb_trial->keep_stats_type( moves::all_stats );
	rb_trial->keep_stats_type( accept_reject );

	RepeatMoverOP rb_cycle( new RepeatMover( rb_trial, inner_cycles_ ) );
	rb_cycle->apply( pose );

	pose = mc_->lowest_score_pose();
	mc_->reset( pose );

	accept_rate_ = rb_trial->acceptance_rate();
}

moves::MonteCarloOP DockingLowRes::get_mc() { return mc_; }


/// @details  Show the complete setup of the docking protocol
void
DockingLowRes::show( std::ostream & out ) const {
	using namespace ObjexxFCL::format;

	// All output will be 80 characters - 80 is a nice number, don't you think?
	std::string line_marker = "///";
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	out << line_marker << A( 47, " Docking Low Res Protocol" ) << space( 27 ) << line_marker << std::endl;
	out << line_marker << space( 74 ) << line_marker << std::endl;

	// Display the number of inner cycles during low res docking
	out << line_marker << " Centroid Inner Cycles: " << inner_cycles_ ;
	out << space( 48 ) << line_marker << std::endl;

	// Display the number of outer cycles during low res docking
	out << line_marker << " Centroid Outer Cycles: " << outer_cycles_;
	out << space( 48 ) << line_marker << std::endl;

	// Display the state of the filters (on or off)
	out << line_marker << " Scorefunction: " << space( 58 ) << line_marker << std::endl;
	scorefxn_->show(out);
	out <<std::endl;

	// Close the box I have drawn
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
}

std::ostream & operator<<(std::ostream& out, const DockingLowRes & dp)
{
	dp.show( out );
	return out;
}

} // namespace docking
} // namespace protocols
