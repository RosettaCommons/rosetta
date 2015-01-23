// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/abinitio/abscript/AbscriptStageMover.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/abinitio/abscript/AbscriptStageMover.hh>

// Package headers
#include <protocols/abinitio/abscript/StagePreparer.hh>

// Project headers
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/pose/Pose.hh>

#include <protocols/moves/MoverFactory.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>

// Utility Headers
#include <utility/vector0.hh>

// Basic Headers
#include <basic/prof.hh>
#include <basic/Tracer.hh>

//option includes
#include <basic/options/option.hh>
#include <basic/options/keys/fold_cst.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>

// C++ Headers
#include <algorithm>

// ObjexxFCL Headers

static thread_local basic::Tracer tr( "protocols.abinitio.abscript.AbscriptStageMover", basic::t_info );

namespace protocols {
namespace abinitio {
namespace abscript{

using namespace core::environment;

AbscriptStageMover::AbscriptStageMover( StageID stage,
                                        moves::MonteCarloOP mc,
                                        core::scoring::ScoreFunctionOP score,
                                        core::Size cycles ):
  stage_( stage ),
  cycles_( cycles ),
  score_( score ),
  mc_( mc ),
  temperature_( 2.0 ),
  seqsep_slope_( 0.0 ),
  seqsep_intcpt_( 0.0 ),
  chbreak_slope_( 0.0 ),
  chbreak_intcpt_( 0.0 )
{
  random_mover_ = moves::RandomMoverOP( new moves::RandomMover() );
}

std::string AbscriptStageMover::get_name() const {
  return "AbscriptStageMover"+utility::to_string( stage_ );
}

void AbscriptStageMover::yield_submovers( MoverSet& set ) const {
  std::copy( preparers_.begin(), preparers_.end(), std::inserter( set, set.begin() ) );
  std::copy( movers_.begin(), movers_.end(), std::inserter( set, set.begin() ) );
}

void AbscriptStageMover::add_submover( ClientMoverOP mover, core::Real weight ) {
  random_mover_->add_mover( mover, weight );
  movers_.insert( mover );
}

void AbscriptStageMover::update_max_seq_sep( core::pose::Pose& pose, core::Real const& progress ){
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  assert( progress >= 0.0 && progress <= 1.0 );

  if( option[ jumps::ramp_chainbreaks ]){
    core::Real const increase_factor = option[ jumps::increase_chainbreak ];
    core::Real new_weight = ( chbreak_slope_ * progress ) + chbreak_intcpt_;

    score_->set_weight( core::scoring::linear_chainbreak, new_weight * increase_factor );
    if ( option[ fold_cst::ramp_coord_cst ] ) {
      score_->set_weight( core::scoring::coordinate_constraint, new_weight * increase_factor );
    }
    if ( option[ jumps::overlap_chainbreak ] && ( stage_ == IVa || stage_ == IVb ) ) {
      score_->set_weight( core::scoring::overlap_chainbreak, progress * increase_factor );
    }

    tr.Info << "chainbreak scores set to " << new_weight << " and ";
  }

  // Configure sequence separation constraints for this run
  constraints_ = constraints_additional::MaxSeqSepConstraintSetOP( new constraints_additional::MaxSeqSepConstraintSet( *pose.constraint_set(), pose.fold_tree()) );
  pose.constraint_set( constraints_ );

  //std::min nullifies 1.2 in setting, but this is how it is in the old code...
  core::Real new_seq_sep_factor = std::min( ( seqsep_slope_ * progress ) + seqsep_intcpt_,
                                            1.0 );

  core::Size max_ft_dist = constraints_->largest_possible_sequence_sep( pose );
  assert( max_ft_dist <= pose.total_residue() );

  tr.Debug << "using factor " << new_seq_sep_factor << " with max ft dist " << max_ft_dist << std::endl;
  core::Size seqsep = new_seq_sep_factor * max_ft_dist;
  tr.Info << "maximum constraint/chainbreak 'seqsep' set to " << seqsep << "." << std::endl;

  if( constraints_->max_seq_sep() != seqsep ){
    constraints_->set_max_seq_sep( seqsep );
  }

  mc_->recover_low( pose );

  // Replace the score function, update maximum sequence sparation for chainbreak scores
  core::scoring::methods::EnergyMethodOptions updated_options( score_->energy_method_options() );
  updated_options.cst_max_seq_sep( seqsep ); //sets both LinearChainrbreak and constraint seqsep
  core::scoring::ScoreFunctionOP new_scorefxn = score_->clone();
  new_scorefxn->set_energy_method_options(updated_options);
  score_ = new_scorefxn;
}

core::Real score_cst(core::pose::Pose& pose, core::scoring::ScoreFunction const& scfxn) {
  scfxn( pose ); //accumulate energies
  return pose.energies().total_energies()[ core::scoring::atom_pair_constraint  ] +
         pose.energies().total_energies()[ core::scoring::coordinate_constraint ] +
         pose.energies().total_energies()[ core::scoring::dihedral_constraint   ] +
         pose.energies().total_energies()[ core::scoring::angle_constraint      ];
}

bool AbscriptStageMover::setup_stage( core::pose::Pose& pose, core::Real const& progress ){

  core::Real old_cst_score = score_cst( pose, *score_ );
  update_max_seq_sep( pose, progress );

  // Skip this stage if changing max_seq_sep doesn't do anything.
  if( stage_ == I && progress > 0.0 &&
     std::abs( old_cst_score - score_cst( pose, *score_ ) ) < 0.01 ){
    return false;
  }

  //Copies the set-up score function into the mc object.
  mc_->score_function( *score_ );

  if( stage_ == I ){
    mc_->set_autotemp( false, temperature_ );
  } else {
    mc_->set_autotemp( true, temperature_ );
  }

  mc_->reset( pose );

  for( PreparerSet::const_iterator preparer = preparers_.begin();
      preparer != preparers_.end(); ++preparer ){
    (*preparer)->prepare( pose, progress );
  }

  ( *score_ )( pose );
  return true;
}

void AbscriptStageMover::apply( core::pose::Pose& pose ){
  using namespace moves;
  using namespace core::scoring::constraints;

  ConstraintSetOP orig_constraints( pose.constraint_set()->clone() );
  mc_->reset( pose );

  RepeatMover( MoverOP( new TrialMover( random_mover_, mc_ ) ), cycles_*cycles_adjust_ ).apply( pose );

  if( stage_ == I ){
    mc_->reset( pose );
  }

  mc_->recover_low( pose );

  if( orig_constraints ){
    pose.constraint_set( orig_constraints );
    ( *score_ )( pose ); //adding constraints zeroes the energies.
  }
}

void AbscriptStageMover::add_preparer( StagePreparerOP mover ) {
  preparers_.insert( mover );
}

void AbscriptStageMover::set_cycles_adjust( core::Real in ){
  cycles_adjust_ = in;
}

void AbscriptStageMover::set_seq_sep_ramping( core::Real const slope, core::Real const intercept ){
  seqsep_slope_ = slope;
  seqsep_intcpt_ = intercept;
}

void AbscriptStageMover::set_chainbreak_ramping( core::Real const slope, core::Real const intercept ){
  chbreak_slope_ = slope;
  chbreak_intcpt_ = intercept;
}



} // abscript
} // abinitio
} // protocols
