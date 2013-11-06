// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (C) 199x-2008 Hebrew University, Jerusalem
//
/// @file   FlexPepDockingLowRes.hh
///
/// @brief low-resolution part of docking protocol
/// @date August 5, 2008
/// @author Barak Raveh


#include <protocols/flexpep_docking/FlexPepDockingFlags.hh>
#include <protocols/flexpep_docking/FlexPepDockingLowRes.hh>

#include <core/types.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <core/fragment/FragSet.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <utility/exit.hh>
#include <string>
// AUTO-REMOVED #include <cstdio>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/chemical/ChemicalManager.fwd.hh>

using namespace protocols::flexpep_docking;

static basic::Tracer TR("protocols.flexPepDockingLowRes");


//////////////////////////////////////////////
/// @brief
/// constructor for low resolution flexpible peptide docking
//
// @param[in] scorefxn_in
//            The scoring function used for optimization
// @param[in] rb_jump
//            The FoldTree rigid body jump over
//            which rigid-body pertrubations are made
FlexPepDockingLowRes::FlexPepDockingLowRes
( FlexPepDockingFlags flags_in,
  core::scoring::ScoreFunctionOP scorefxn_in,
  core::kinematics::MoveMapOP movemap_in,
  Size const rb_jump_in )
  : flags_(flags_in),
    movemap_(movemap_in),
    rb_jump_(rb_jump_in)
{
  using namespace basic::options;

  // TODO: create a set_defaults() function
  scorefxn_ = scorefxn_in->clone();

  // Loop modeling options
  // NOTE: most LoopRelax options are initiated automatically from cmd-line
  // TODO: LoopRelaxMover is a wrapper, perhaps user the LoopModel class explicitly
  loop_relax_mover_ = new protocols::comparative_modeling::LoopRelaxMover();
  // loop_relax_mover_->centroid_scorefxn(scorefxn_); // TODO: we need a chain brteak score here, so let's leave it for modeller default?
  loop_relax_mover_->refine("no"); // centroid modeling only
  loop_relax_mover_->relax("no"); // centroid modeling only
  if( option[ OptionKeys::loops::frag_files ].user() )
    {
      // these protocols optionally take a fragment set .. only load if
      // specified
      utility::vector1< core::fragment::FragSetOP > frag_libs;
      protocols::loops::read_loop_fragments( frag_libs );
      loop_relax_mover_->frag_libs( frag_libs );
    }
}


// empty destructor - for good inclusion of OP clasesses
FlexPepDockingLowRes::~FlexPepDockingLowRes()
{}


///////////////////////////////////////////////
/// @brief initial setup for apply
void
FlexPepDockingLowRes::setup_for_apply( core::pose::Pose& pose )
{
  double temperature = 0.8;
  mc_ = new moves::MonteCarlo( pose, *scorefxn_, temperature );
  // setup minimizer
  std::string min_type = "dfpmin_atol"; // armijo_nonmonotone? different tolerance?
  double min_func_tol = 0.1;
  minimizer_ = new protocols::simple_moves::MinMover(
    movemap_, scorefxn_, min_type, min_func_tol, true /*nb_list accel.*/ );
}


///////////////////////////////////////////////
// switch pose to centroid mode, if not already there
void
FlexPepDockingLowRes::to_centroid
( core::pose::Pose & pose ) const
{
  if(!pose.is_fullatom())
    return;
  TR.Debug << "Switching to centroid" << std::endl;
  protocols::simple_moves::SwitchResidueTypeSetMover
    to_centroid_mover( core::chemical::CENTROID );
  to_centroid_mover.apply(pose);
}


//////////////////////////////////////////////
// switch pose to full-atom mode, if not already
// using side-chains of referencePose
void
FlexPepDockingLowRes::to_allatom
( core::pose::Pose & pose, core::pose::Pose& referencePose ) const
{
  runtime_assert(referencePose.is_fullatom());
  protocols::simple_moves::SwitchResidueTypeSetMover
    to_all_atom_mover( core::chemical::FA_STANDARD );
  protocols::simple_moves::ReturnSidechainMover
    recover_sidechains( referencePose );
  recover_sidechains.apply( pose );
}


////////////////////////////////////////////
void
FlexPepDockingLowRes::torsions_monte_carlo
( core::pose::Pose & pose,
  const int cycles,
  double& acceptance_rate )
{
  using namespace protocols::moves;
  // setup sub-moves
  simple_moves::SmallMoverOP small_mover =
    new protocols::simple_moves::SmallMover( movemap_, 0.8 /*temp*/, 100 /*nmoves ???*/ );
  small_mover->angle_max('L',flags_.smove_angle_range);
  protocols::simple_moves::ShearMoverOP shear_mover =
    new protocols::simple_moves::ShearMover( movemap_, 0.8 /*temp*/, 100 /*nmoves ???*/ );
  shear_mover->angle_max('L',flags_.smove_angle_range);

  // setup cycle of sub-moves
  CycleMoverOP cyclic_tor_moves =
    new CycleMover;
  cyclic_tor_moves->add_mover(small_mover);
  cyclic_tor_moves->add_mover(shear_mover);

  // wrap with monte-carlo trial mover
  TrialMoverOP mc_trial = new TrialMover( cyclic_tor_moves, mc_ );
  mc_trial->keep_stats_type( accept_reject ); // track stats (for acceptance rate)

  // Do initial forced perturbation // TODO: do we want to keep this part?
  small_mover->apply(pose);
  shear_mover->apply(pose);
  // run Monte-Carlo
  for( int i=1; i<=cycles; ++i ) {
    mc_trial->apply( pose );
  }
  // extract best pose and return statistics
  pose = mc_->lowest_score_pose();
  mc_->reset( pose );
  acceptance_rate = mc_trial->acceptance_rate();
}


///////////////////////////////////////////
void
FlexPepDockingLowRes::loopclosure_monte_carlo
( core::pose::Pose & pose )
{
  using namespace protocols::moves;

  if(flags_.peptide_nres() < 5) // TODO: 3 is minimum for KIC loop closure, and flanks should be excluded
    return;
  // set up and model a random loop
  Size first_res = flags_.peptide_first_res() + 1;
  Size last_res = flags_.peptide_last_res() - 1;
  protocols::loops::LoopsOP loops = new protocols::loops::Loops();
  loops->add_loop(first_res, last_res); // TODO: cut defaults to zero, is this a random cut?
  for(Size i = first_res; i <= last_res ; i++)
    runtime_assert( movemap_->get_bb(i) ); // verify loop is movable, this should have been always true for the peptide
  loop_relax_mover_->loops( loops );
  loop_relax_mover_->apply( pose );

}



///////////////////////////////////////////
// returns the acceptance rate from the RB monte-carlo
void
FlexPepDockingLowRes::rigidbody_monte_carlo
( core::pose::Pose & pose,
  const int cycles,
  const float trans_magnitude,
  const float rot_magnitude,
  double& acceptance_rate
)
{
  using namespace protocols::moves;

  // set up monte-carlo trial moves
  rigid::RigidBodyPerturbNoCenterMoverOP rb_mover =
    new rigid::RigidBodyPerturbNoCenterMover(
      rb_jump_, rot_magnitude, trans_magnitude );
  TrialMoverOP mc_trial = new TrialMover( rb_mover, mc_ );
  mc_trial->keep_stats_type( accept_reject ); // track stats (for acceptance rate)

  // run monte-carlo
  for(int i = 1; i <= cycles ; i++)
    {
      mc_trial->apply(pose);
    }

  // extract best pose and return statistics
  pose = mc_->lowest_score_pose();
  mc_->reset( pose );
  acceptance_rate = mc_trial->acceptance_rate();
}


///////////////////////////////////////////////
void
FlexPepDockingLowRes::apply( core::pose::Pose & pose )
{
  using namespace core;
  // variable declarations: //
  pose::Pose startPose = pose;
  int outer_cycles=10; // TODO: runtime param?
  int inner_cycles_rb = 50; // TODO: runtime param?
  int inner_cycles_torsions = 50; // TODO: runtime param?
  double trans_mag = 0.5; // Angstrom // TODO: runtime param?
  double rot_mag = 5; // Degrees // TODO: runtime param?
  double rb_acceptance, torsions_acceptance; // MC acceptance rates
  // real stuff: //
  to_centroid(pose);
  setup_for_apply(pose);
  TR.flush();
  TR.flush();
  for ( int i=1; i<=outer_cycles; ++i) {
    if(flags_.rbMCM) {
      rigidbody_monte_carlo
	( pose , inner_cycles_rb, trans_mag, rot_mag, rb_acceptance);
      minimizer_->apply(pose);
    }
    if(flags_.torsionsMCM) {
      torsions_monte_carlo
	( pose , inner_cycles_torsions, torsions_acceptance);
      minimizer_->apply(pose);
    }
    if(flags_.peptide_loop_model){
      loopclosure_monte_carlo( pose );
      minimizer_->apply(pose);
    }
  }
  TR.flush();
  to_allatom(pose, startPose /*reference pose*/);
}

std::string
FlexPepDockingLowRes::get_name() const {
	return "FlexPepDockingLowRes";
}

