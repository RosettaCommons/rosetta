
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @detailed
/// @author Oliver Lange
///


// Unit Headers
#include <protocols/loops/loop_closure/ccd/WidthFirstSlidingWindowLoopClosure.hh>

// Package Headers
// AUTO-REMOVED #include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/loop_closure/ccd/LoopClosure.hh>
#include <protocols/loops/loop_closure/ccd/ShortLoopClosure.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Exceptions.hh>

#include <protocols/evaluation/util.hh>
// Project Headers
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/util.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>

#include <core/fragment/FrameList.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>
// AUTO-REMOVED #include <core/fragment/FragCache.hh> // for FragStore

// AUTO-REMOVED

#include <basic/options/option.hh> // for quick-test from run:dry_run
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>
//#include <basic/options/keys/out.OptionKeys.gen.hh> //for frag_file ... temporary!
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/fast_loops.OptionKeys.gen.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>

#include <core/fragment/FragData.hh>
#include <utility/vector1.hh>

//Auto Headers


//numeric headers


static thread_local basic::Tracer tr( "protocols.loops.loop_closure.ccd.WidthFirstSlidingWindowLoopClosure" );

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

using namespace core;
using namespace pose;
using namespace fragment;

WidthFirstSlidingWindowLoopClosure::WidthFirstSlidingWindowLoopClosure(
  fragment::FragSetCOP fragset,
  scoring::ScoreFunctionOP scorefxn,
  kinematics::MoveMapCOP movemap
) : SlidingWindowLoopClosure( fragset, scorefxn, movemap )
{
  set_defaults();
}

WidthFirstSlidingWindowLoopClosure::WidthFirstSlidingWindowLoopClosure() :
  SlidingWindowLoopClosure()
{
  set_defaults();
}

std::string
WidthFirstSlidingWindowLoopClosure::get_name() const {
	return "WidthFirstSlidingWindowLoopClosure";
}

void WidthFirstSlidingWindowLoopClosure::set_defaults() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  Parent::set_defaults();
  window_acceptance_ratio_ = option[ fast_loops::window_accept_ratio ]();
  nr_scored_sampling_passes_ = option[ fast_loops::nr_scored_sampling_passes ]();
  min_breakout_good_loops_ = option[ fast_loops::min_breakout_good_loops ]();
  min_breakout_fast_loops_ = option[ fast_loops::min_breakout_fast_loops ]();
  min_good_loops_ = option[ fast_loops::min_good_loops ]();
  min_fast_loops_ = option[ fast_loops::min_fast_loops ]();
  vdw_delta_ = option[ fast_loops::vdw_delta ]();
  give_up_ = option[ fast_loops::give_up ]();
  nr_scored_fragments_ = option[ fast_loops::nr_scored_fragments ]();
  chainbreak_max_ = basic::options::option[ basic::options::OptionKeys::fast_loops::chainbreak_max ]();
  tr.Info << "WidthFirstSlidingWindowLoopClosure::defaults " << std::endl;
}

void WidthFirstSlidingWindowLoopClosure::register_options() {
    using namespace basic::options;
	using namespace basic::options::OptionKeys;

    option.add_relevant(  fast_loops::window_accept_ratio );
    option.add_relevant(  fast_loops::nr_scored_sampling_passes );
    option.add_relevant(  fast_loops::nr_scored_fragments );
    option.add_relevant(  fast_loops::min_breakout_good_loops );
    option.add_relevant(  fast_loops::min_breakout_fast_loops );
    option.add_relevant(  fast_loops::min_good_loops );
    option.add_relevant(  fast_loops::min_fast_loops );
    option.add_relevant(  fast_loops::vdw_delta );
    option.add_relevant(  fast_loops::give_up );
}

void
WidthFirstSlidingWindowLoopClosure::sample_loops( Pose& more_cut, Pose& less_cut ) {
  if ( bQuickTest() ) return; //let's say we found a good loop

  if( basic::options::option[ basic::options::OptionKeys::run::test_cycles ]() ){
    min_loop_size_ = 6;
    max_loop_size_ = 7;
    min_good_loops_ = 1;
    scored_frag_cycle_ratio_ = 0.05;
    short_frag_cycle_ratio_ = 0.05;
    vdw_delta_ = 200;
    score_delta_ = 200;
  }

  const Real REALLY_BAD_SCORE ( 1000000000.0 );
  best_score_ = REALLY_BAD_SCORE;

  if ( bKeepFragments_ ) {
    closure_fragments_ = new fragment::OrderedFragSet;
  }


  tr.Debug << "Trying loop-sizes: " << loop_ << std::endl;
  Size const actual_max_loop_size( loop_.size() );
  Size min_breakout_good_loops( std::max( min_breakout_good_loops_, min_good_loops_ ) );
  tr.Info << "breakout good loops " << min_breakout_good_loops_ << " min_good_loops_ " << min_good_loops_
	  << "breakout fast loops " << min_breakout_fast_loops_ << " min_fast_loops_ " << min_fast_loops_ << std::endl;
  tr.Debug << "Trying loop-sizes: "
	   << " "   << std::min( actual_max_loop_size, min_loop_size_ )
	   << " " 	 << " - "
	   << " "	 << std::min( actual_max_loop_size, max_loop_size_ )
	   << " " 	 << std::min( actual_max_loop_size, max_loop_size_ ) -  std::min( actual_max_loop_size, min_loop_size_ )
	   << " "   << std::endl;

  tr.Debug << "LOOP: " << loop_ << std::endl;

  Size fast_loop_count( 0 );
  Size good_loop_count( 0 );
  Size scored_frags_attempted( 0 );

  //we need to move this into Base class ... right now it has effect to set scorefxn_ as well...
  scoring::ScoreFunctionOP frag_scorefxn = setup_frag_scorefxn();
  tr.Debug << "Trying loop-sizes: " << loop_ << std::endl;
  tr.Info << "---------------- LOOP SAMPLING based on this scorefunction: ----------------\n";
  if ( tr.Info.visible() ) frag_scorefxn->show( tr.Info, more_cut );
  tr.Info << std::endl;

  tr.Debug << "Trying loop-sizes: " << loop_ << std::endl;
  tr.Info << "---------------- LOOP SELECTION based on this scorefunction: ----------------\n";
  if ( tr.Info.visible() ) scorefxn_->show( tr.Info, more_cut );
  tr.Info << std::endl;

  loops::remove_cutpoint_variants( more_cut, true );
  loops::add_single_cutpoint_variant( more_cut, loop_ );

  loops::remove_cutpoint_variants( less_cut, true  );
  loops::add_single_cutpoint_variant( less_cut, loop_ );

  tr.Debug << "MOREFOLDTREE: " << more_cut.fold_tree();
  tr.Debug << "LESSFOLDTREE: " << less_cut.fold_tree();
  if ( evaluator_ && tr.Debug.visible() ) evaluate_pose( more_cut, *evaluator_, tr.Debug );
  if ( evaluator_ && tr.Debug.visible() ) evaluate_pose( less_cut, *evaluator_, tr.Debug );

  {	//take initial loop-conformation as closing candidate
    tr.Debug << "CAPTURE INITIAL POSE LOOPS" << std::endl;
    using namespace fragment;
    FrameOP closure_frame = new Frame( loop_.start(), new FragData( new BBTorsionSRFD, loop_.size() ) );
    FrameList closure_frames;
    closure_frame->steal( more_cut );

    CCDLoopClosureMover fast_ccd( loop_, movemap_ );
    fast_ccd.apply( more_cut );
    if ( fast_ccd.success() ) {
      closure_frame->steal( more_cut );
    }

    closure_frames.push_back( closure_frame );

    good_loop_count+=process_fragments( closure_frames, more_cut, less_cut );
    fast_loop_count = good_loop_count;
    tr.Debug << "INITIAL POSE yielded " << good_loop_count << " valuable conformations " << std::endl;
  } //scope

  //	pose::Pose best_pose = more_cut;
    fragment::FrameList closure_frames;

  // try different loop sizes
  WindowList windows;
  for ( Size loop_size = std::min( actual_max_loop_size, min_loop_size_ );
	loop_size <= std::min( actual_max_loop_size, max_loop_size_ ); ++loop_size ) {
    tr.Debug << "loop-size: " << loop_size << std::endl;
    // try different sliding windows, sorted by loop_fraction in window
    generate_window_list( loop_size, windows ); //generates a list of windows sorted by loop-fraction
  }
  windows.sort();
  windows.reverse();

  WindowList probable_windows; //fill these with windows that give success with the fast ShortLoopClosure


  for ( WindowList::const_iterator it = windows.begin(), eit = windows.end();
	it != eit && fast_loop_count < min_breakout_fast_loops_ && good_loop_count < min_breakout_good_loops;
	++it ) {

    fragment::FrameList closure_frames;
    // close this loop
    Loop const current_loop = it->second;
    tr.Debug << "attempt closure on " << current_loop << std::endl;

    tr.Info << " unscored (vdw only) fragment sampling... " << current_loop <<  std::endl;
    ShortLoopClosure fast_closure( fragset_, current_loop, movemap_ );
    fast_closure.set_cycles( short_frag_cycle_ratio_ ); //determines number of closure attempts and frag_cycles per attempt
    fast_closure.ramp_chainbreak();
    if ( fast_closure.apply( more_cut ) ) {
      closure_frames.push_back( fast_closure.closure_fragments() );
      Size new_loops = process_fragments( closure_frames, more_cut, less_cut );
      if ( 1.0*new_loops > ( window_acceptance_ratio_ * fast_closure.nr_fragments() ) ) {
	tr.Info << "good window keep this for scored sampling stage" << std::endl;
	probable_windows.push_back( *it );
      }
      fast_loop_count += new_loops;
    }
  }

  //I don't really see why this is useful... also with min_fast_loops_ = 0 this could be problematic...
  //if ( fast_loop_count < min_fast_loops_ ) {
  //  min_breakout_good_loops = min_fast_loops_; //if we don't have enough fast loops -- try to get at least as many good loops...
  //}

  tr.Info << "finished unscored sampling: loops found: " << fast_loop_count << " (needs " << min_fast_loops_
	  << " ) good windows: " << probable_windows.size() << " min_breakout_good_loops " << min_breakout_good_loops << std::endl;
  for ( Size attempt_count = 0; attempt_count < nr_scored_sampling_passes_; attempt_count++ ) {
    if ( good_loop_count >= min_breakout_good_loops ) break;
    if ( scored_frags_attempted > give_up_ && good_loop_count == 0 ) break;

    for ( WindowList::const_iterator it = probable_windows.begin(),
	    eit = probable_windows.end(); it != eit; ++it ) {

      if ( good_loop_count >= min_breakout_good_loops ) break;
      if ( scored_frags_attempted > give_up_ && good_loop_count == 0 ) break;
      scored_frags_attempted += nr_scored_fragments_;

      fragment::FrameList closure_frames;

      Loop const current_loop = it->second;

      tr.Info << "scored fragment sampling on ... " << current_loop <<  std::endl;
      LoopClosure scored_closure( fragset_, frag_scorefxn, current_loop, movemap_ );
      scored_closure.set_cycles( scored_frag_cycle_ratio_ );
      scored_closure.set_nr_fragments( nr_scored_fragments_ ); //make only 5 fragment per window ... but reiterate windows a couple of times
      scored_closure.ramp_chainbreak();
      //		if ( bIdealLoopClosing() ){
      // note this apply doesn't change the more_cut pose

      if ( scored_closure.apply( more_cut ) ) closure_frames.push_back( scored_closure.closure_fragments() );

      good_loop_count += process_fragments( closure_frames, more_cut, less_cut );
      tr.Info << "process fragments... GoodLoops: " <<  good_loop_count << std::endl;

    } // probable_windows    //		if ( bIdealLoopClosing() ){
  } // for ( attempt_count ) ... window-scans

  tr.Info << "finished sampling... Good scored Loops: " <<  good_loop_count << " required: " << min_good_loops_
	  << "good loops: " << fast_loop_count+good_loop_count << " ( " << min_fast_loops_ << " ) " <<  std::endl;

  if ( !best_fragment_.is_valid() || ( good_loop_count < min_good_loops_ && (fast_loop_count + good_loop_count) < min_fast_loops_ ) ) {
    tr.Warning << "WARNING: no good loop found !" << std::endl;
    throw( EXCN_Loop_not_closed() );
  }
} //sample_loops

} // namespace ccd
} // namesapce loop_closure
} // namespace loops
} // namespace protocols
