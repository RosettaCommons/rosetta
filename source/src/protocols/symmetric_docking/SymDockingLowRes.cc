// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DockingLowRes
/// @brief protocols that are specific to docking low resolution
/// @details This is to a very large extent a copy of the docking
/// @details protocol. Should derive out of that class instead.
/// @author Ingemar Andre

#include <protocols/symmetric_docking/SymDockingLowRes.hh>

// Rosetta Headers
#include <core/kinematics/MoveMap.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/OutputMovers.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/pose/symmetry/util.hh>
#include <core/kinematics/FoldTree.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <string>

//Utility Headers

#include <basic/options/keys/docking.OptionKeys.gen.hh>

#include <numeric/trig.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>

#include <basic/Tracer.hh>

#include <protocols/moves/MoverContainer.hh>
#include <utility>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.symetric_docking.SymDockingLowRes" );

//     originally from dock_structure.cc Jeff Gray April 2001
//

using namespace core;

namespace protocols {
namespace symmetric_docking {

// constructor with arguments
SymDockingLowRes::SymDockingLowRes(
	core::scoring::ScoreFunctionCOP scorefxn_in
) : Mover(), scorefxn_(std::move(scorefxn_in))
{
	moves::Mover::type( "SymDockingLowRes" );
}

SymDockingLowRes::~SymDockingLowRes()= default;

moves::MoverOP
SymDockingLowRes::clone() const {
	return moves::MoverOP( new SymDockingLowRes( *this ) );
}

void
SymDockingLowRes::set_default( core::pose::Pose & pose ) {
	using namespace basic::options;

	// sets up the stuff in pose
	(*scorefxn_)( pose );

	// cycles
	inner_cycles_ = option[ OptionKeys::docking::docking_centroid_inner_cycles ]();
	outer_cycles_ = option[ OptionKeys::docking::docking_centroid_outer_cycles ]();

	if ( option[ OptionKeys::docking::dock_mcm_trans_magnitude ].user() ) {
		trans_magnitude_ = option[ OptionKeys::docking::dock_mcm_trans_magnitude ]();
	} else {
		trans_magnitude_ = 1.5;
	}

	if ( option[ OptionKeys::docking::dock_mcm_rot_magnitude ].user() ) {
		rot_magnitude_ = option[ OptionKeys::docking::dock_mcm_rot_magnitude ]();
	} else {
		rot_magnitude_ = 4;
	}

	chi_ = false;
	bb_ = false;

	temperature_ = 0.8;

	nb_list_ = true; /// not sure if this should be true or not
	accept_rate_ = 0.0;

	set_default_mc( pose );
	set_default_move_map( pose );
	set_default_protocol( pose );
}

moves::MonteCarloOP
SymDockingLowRes::get_mc() { return mc_; }

void
SymDockingLowRes::set_default_mc( pose::Pose & pose ) {
	// create the monte carlo object and movemap
	mc_ = moves::MonteCarloOP( new moves::MonteCarlo( pose, *scorefxn_, temperature_ ) );
}

void SymDockingLowRes::set_default_move_map( pose::Pose & pose ) {
	using namespace core::conformation::symmetry;

	movemap_ = core::kinematics::MoveMapOP( new kinematics::MoveMap() );
	movemap_->set_bb( bb_ );
	movemap_->set_chi( chi_ );
	core::pose::symmetry::make_symmetric_movemap( pose, *movemap_ );

}

void SymDockingLowRes::set_default_protocol( pose::Pose & pose ){
	using namespace moves;
	using namespace conformation::symmetry;

	debug_assert( core::pose::symmetry::is_symmetric( pose ));
	auto & symm_conf (
		dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

	std::map< Size, SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );

	rb_mover_ = rigid::RigidBodyDofSeqPerturbMoverOP( new rigid::RigidBodyDofSeqPerturbMover( dofs , rot_magnitude_, trans_magnitude_ ) );

	docking_lowres_protocol_ = moves::SequenceMoverOP( new SequenceMover );
	docking_lowres_protocol_->add_mover( rb_mover_ );

	if ( basic::options::option[basic::options::OptionKeys::docking::multibody].user() ) {
		utility::vector1<int> mbjumps = basic::options::option[basic::options::OptionKeys::docking::multibody]();
		for ( Size ij = 1; ij <= symm_conf.Symmetry_Info()->get_njumps_subunit(); ++ij ) {
			if ( mbjumps.size()==0 || std::find(mbjumps.begin(),mbjumps.end(),(int)ij)!=mbjumps.end() ) {
				TR << "add subunit jump mover " << ij << std::endl;
				docking_lowres_protocol_->add_mover( MoverOP( new rigid::RigidBodyPerturbMover(ij,rot_magnitude_,trans_magnitude_) ) );
			}
		}
	}

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
void SymDockingLowRes::apply( core::pose::Pose & pose )
{
	using namespace scoring;

	TR << "in DockingLowRes.apply\n";

	set_default( pose );

	TR << "::::::::::::::::::Centroid Rigid Body Adaptive:::::::::::::::::::\n";

	for ( int i=1; i<=outer_cycles_; ++i ) {
		rigid_body_trial( pose );
		if ( accept_rate_ < 0.5 ) {
			trans_magnitude_ *= 0.9;
			rot_magnitude_ *= 0.9;
		} else {
			trans_magnitude_ *= 1.1;
			rot_magnitude_ *= 1.1;
		}
		// if ( jump_out_check() ) return;
	}
	mc_->recover_low( pose );
	TR.flush();
	//pose.energies().show( std::cout );
}

std::string
SymDockingLowRes::get_name() const {
	return "SymDockingLowRes";
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
void SymDockingLowRes::rigid_body_trial( core::pose::Pose & pose )
{
	using namespace moves;

	PDBDumpMoverOP dump( new PDBDumpMover("lowres_cycle_") );
	// dump->apply( pose );
	MCShowMoverOP mc_show( new MCShowMover( mc_ ) );
	// mc_show->apply( pose );

	rb_mover_->rot_magnitude( rot_magnitude_ );
	rb_mover_->trans_magnitude( trans_magnitude_ );

	TrialMoverOP rb_trial( new TrialMover( docking_lowres_protocol_, mc_ ) );
	rb_trial->keep_stats_type( moves::accept_reject );

	RepeatMoverOP rb_cycle( new RepeatMover( rb_trial, inner_cycles_ ) );

	rb_cycle->apply( pose );

	pose = mc_->lowest_score_pose();
	//pose.energies().show( std::cout );
	mc_->reset( pose );

	accept_rate_ = rb_trial->acceptance_rate();
}

} // namespace docking
} // namespace protocols
