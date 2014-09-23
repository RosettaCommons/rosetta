// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PoseEvaluator
/// @brief PoseEvaluator
/// @detailed
///
///
/// @author Oliver Lange



// Unit Headers
#include <protocols/simple_filters/EvaluatedTrialMover.hh>

// Package Headers
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/evaluation/PoseEvaluator.hh>

// Project Headers
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>

// ObjexxFCL Headers

// Utility headers
// AUTO-REMOVED #include <basic/Tracer.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/prof.hh>

#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/simple_filters/ScoreEvaluator.hh>
#include <utility/vector1.hh>

// C++ headers


namespace protocols {
namespace simple_filters {

using namespace core;
using namespace core::io::silent;

/// c'stor
EvaluatedTrialMover::EvaluatedTrialMover(
	moves::MoverOP mover_in,
	moves::MonteCarloOP mc_in,
	evaluation::PoseEvaluatorOP evaluator_in,
	std::string tag
) :
	TrialMover( mover_in, mc_in),
	tag_( tag )
{
	using protocols::evaluation::PoseEvaluatorOP;
	evaluator_ = evaluation::MetaPoseEvaluatorOP( new evaluation::MetaPoseEvaluator );
	evaluator_->add_evaluation( PoseEvaluatorOP( new simple_filters::ScoreEvaluator( "full",mc_in->score_function().clone() ) ) );
	evaluator_->add_evaluation( evaluator_in );
}

EvaluatedTrialMover::~EvaluatedTrialMover() {}

void
EvaluatedTrialMover::apply( pose::Pose& pose ) {
	/// OKAY this IS code duplication but it is just a HACK
	/// it's true. code duplication is totally ok when it's just a hack. :P
	using namespace protocols::moves;

  using scoring::total_score;
  /// get the initial scores
  if ( keep_stats_type() == all_stats ) {
    stats_.add_score( mc_->last_accepted_score() ); ///< initial_last_accepted_score
    stats_.add_score( pose.energies().total_energy() ); ///< initial_pose_score
  }

  /// make the move
  mover_->apply( pose );

  PROF_START( basic::TEST3 );
  SilentStructOP pss( new io::silent::ProteinSilentStruct );
  evaluator_->apply( pose, tag_, *pss );
  PROF_STOP( basic::TEST3 );

  if ( keep_stats_type() == all_stats ) { //// score and get residue energies
    //Real newscore =
    mc_->score_function()( pose );
    //std::cout << "newscore: " << newscore << std::endl;
    /// Now handled automatically.  mc_->score_function().accumulate_residue_total_energies( pose );
    stats_.add_score( pose.energies().total_energy() ); ///< score_after_move
    //if ( newscore != pose.energies().total_energy() ) {
    //	std::cout << " newscore != perceived newscore: " << newscore << " != " <<
    //		pose.energies().total_energy() << std::endl;
    //}
  }

  /// test if MC accepts or rejects it
  bool accepted_move = mc_->boltzmann( pose, mover_->type() );

  PROF_START( basic::TEST3 );
  pss->add_energy( "acceptance", accepted_move );
  evals_.push_back( pss );
  PROF_STOP( basic::TEST3 );

  if ( keep_stats_type() > no_stats ) {
    stats_.accepted( accepted_move );
    stats_.print( mc_, mover_->type() );
    // std::cout << "Acceptance rate: " << stats_.acceptance_rate() << std::endl;
  }


}

std::string
EvaluatedTrialMover::get_name() const {
	return "EvaluatedTrialMover";
}

void
EvaluatedTrialMover::dump_file( std::string file ) {
  //copy some of the functionality of SilentFileData::write_silent_struct
  // want to avoid opening/closing of file for each line of data ... write as burst
  // might add this to SilentFileData
  if ( evals_.size() ) {
    utility::io::ozstream output;
    if ( !utility::file::file_exists( file ) ) {
      output.open( file );
      evals_[ 1 ]->print_header( output );
    } else {
      output.open_append( file );
    }

    for ( SilentInfoList::const_iterator it=evals_.begin(), eit=evals_.end(); it!=eit; ++it ) {
      (*it)->print_scores( output );
    }
  }
}



}
}
