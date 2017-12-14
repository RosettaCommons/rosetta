// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @details
/// @author Oliver Lange


// Unit Headers
#include <protocols/loops/loop_closure/ccd/LoopClosure.hh>

// Package Headers
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/loops_main.hh>

// Project Headers
#include <core/pose/Pose.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/fragment/Frame.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>


#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialMover.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>

//numeric headers

//// C++ headers
#include <cstdlib>
#include <string>

#include <core/fragment/FragData.hh>
#include <utility>
#include <utility/vector1.hh>


namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

using namespace core;
using namespace pose;

static basic::Tracer tr( "protocols.loops.loop_closure.ccd.LoopClosure" );

LoopClosure::LoopClosure(
	fragment::FragSetCOP fragset,
	scoring::ScoreFunctionOP scorefxn,
	Loop loop_def,
	kinematics::MoveMapCOP movemap
) : loop_ ( loop_def ),
	scorefxn_(std::move( scorefxn )),
	movemap_(std::move( movemap )),
	frag_mover_( /* NULL */ ),
	ccd_mover_( /* NULL */ ),
	fragset_(std::move( fragset )),
	bEnableCcdMoves_( loop_def.size() >= 10 ? true : false ),
	bRampChainbreak_( false )
{
	set_cycles( 1.0 );
	temperature_ = 2.0;
	init();
}

LoopClosure::LoopClosure() :
	scorefxn_( /* NULL */ ),
	movemap_( /* NULL */ ),
	frag_mover_( /* NULL */ ),
	ccd_mover_( /* NULL */ ),
	fragset_( /* NULL */ ),
	bEnableCcdMoves_( false )
{
	set_cycles( 1.0 );
	temperature_ = 2.0;
}

core::kinematics::MoveMapCOP LoopClosure::movemap() const { return movemap_; }

void LoopClosure::set_movemap( core::kinematics::MoveMapCOP mm ) { movemap_ = mm; }

void LoopClosure::set_fragset( core::fragment::FragSetCOP frags ) { fragset_ = frags; }

void LoopClosure::init() {
	runtime_assert( fragset_ != nullptr );
	runtime_assert( movemap_ != nullptr );
	runtime_assert( scorefxn_ != nullptr );

	//make movemap that only allows bb moves within loop ( if master movemap allows the move, too ).
	kinematics::MoveMapOP loop_movemap( new kinematics::MoveMap );
	loop_movemap->set_bb( false );
	for ( Size pos = loop_.start(); pos <= loop_.stop(); pos ++ ) {
		if ( movemap_->get_bb( pos ) ) loop_movemap->set_bb( pos, true );
	}
	set_movemap( loop_movemap );

	simple_moves::ClassicFragmentMoverOP ptr;
	frag_mover_ = ptr = simple_moves::ClassicFragmentMoverOP( new simple_moves::ClassicFragmentMover( fragset_, movemap_ ) );
	ptr->enable_end_bias_check( false ); //uniform sampling
	ptr->set_check_ss( false );
	if ( bEnableCcdMoves_ ) {
		ccd_mover_ = moves::MoverOP( new CCDLoopClosureMover( loop_, movemap_ ) );
	}
	runtime_assert( loop_.size() > 0 );
	init_mc();
	// put a template frame at top of list -- will always catch the fragments according to closure_frames.front()
	using namespace fragment;
	closure_frame_ = core::fragment::FrameOP( new Frame( loop_.start(), FragDataCOP( FragDataOP( new FragData( SingleResidueFragDataOP( new BBTorsionSRFD ), loop_.size() ) ) ) ) );
}

LoopClosure::~LoopClosure() = default;

void LoopClosure::set_cycles( core::Real cycle_ratio ) {
	nr_fragments_ = static_cast< int > (100*cycle_ratio);
	cycles_ = static_cast< int > (20*std::max( (int) loop_.size(), 5 /*arbitrary choice*/ ) *cycle_ratio);
}

void LoopClosure::set_nr_fragments( core::Size nr_fragments ) {
	nr_fragments_ = nr_fragments;
}

void LoopClosure::init_mc() {
	mc_ = moves::MonteCarloOP( new moves::MonteCarlo( *scorefxn_, temperature_ ) );
	final_weight_linear_chainbreak_ = scorefxn_->get_weight( scoring::linear_chainbreak );
	final_weight_overlap_chainbreak_ = scorefxn_->get_weight( scoring::overlap_chainbreak );
}


bool
LoopClosure::apply( pose::Pose const& pose_in ) {
	pose::Pose pose( pose_in );

	loops::add_single_cutpoint_variant( pose, loop() );

	Real best_score( 10000000.0 );
	//Real best_fdev( 100000.0 );
	for ( Size c1 = 1; c1 <= nr_fragments_; ++c1 ) {

		// replace all torsions of loop once
		for ( Size pos = loop_.start(); pos <= (Size) loop_.stop(); pos++ ) {
			bool success ( frag_mover_->apply( pose, pos ) );
			if ( !success ) {
				tr.Warning << " could not make fragment move at loop position " << pos
					<< " seqpos: " << pos << std::endl;
			}
		}

		mc().reset( pose );

		// monte-carlo based minimization of score via the supplied mover ( usually fragment insertions / and evtl. ccd-moves )
		tr.Trace <<" before frag_cycles " << cycles_ << std::endl;
		do_frag_cycles( pose );
		tr.Trace << "after frag_cycles" << std::endl;
		mc().recover_low( pose );
		Real score;
		score = (*scorefxn_)( pose );
		tr.Trace << "start ccd " << std::endl;
		CCDLoopClosureMover fast_ccd( loop_, movemap() );
		fast_ccd.apply( pose );
		Real dev = fast_ccd.deviation();
		//    if ( tr.Trace.visible() ) scorefxn_->show( tr, pose );

		if ( fast_ccd.success() ) {
			best_score = std::min( best_score,score);
			catch_fragment( pose );
		}
		tr.Debug << "LoopClosure: fragment " << c1 << " best_score: " << best_score << " pre ccd score: " <<  score
			<< " dev: " << dev << " " << ( fast_ccd.success() ? "SUCCESS" : "FAIL" )
			<< std::endl;
	}
	return closure_frame_->nr_frags() >= 1;
}

void
LoopClosure::ramp_chainbreak( Size iter, Size total ) const {
	runtime_assert( total > 0 );
	if ( iter < static_cast< Size >( total/2.0) ) {
		core::Real const progress( iter/total * 2.0 );
		scorefxn_->set_weight( scoring::linear_chainbreak, progress * final_weight_linear_chainbreak_ );
		scorefxn_->set_weight( scoring::overlap_chainbreak, 0.0 );
	} else {
		core::Real const progress( iter/total );
		scorefxn_->set_weight( scoring::linear_chainbreak, 1.0 );
		scorefxn_->set_weight( scoring::overlap_chainbreak, progress * final_weight_overlap_chainbreak_ );
	}
	mc_->score_function( *scorefxn_ );
}

void
LoopClosure::do_frag_cycles( pose::Pose &pose ) const {
	moves::TrialMover frag_trial( frag_mover_, mc_ );

	moves::TrialMoverOP ccd_trial;
	if ( ccd_mover_ ) ccd_trial = moves::TrialMoverOP( new moves::TrialMover( ccd_mover_, mc_ ) );
	if ( bRampChainbreak_ ) ramp_chainbreak( Size(1), cycles_ );
	for ( Size i = 1; i <= cycles_; i++ ) {
		frag_trial.apply( pose );
		if ( bRampChainbreak_ && ( i % 20 == 0 ) ) ramp_chainbreak( i, cycles_ );
		if ( i % 10 ==0 ) tr.Trace << "loop-frag-trials: iterations: " << i << std::endl;
		if ( ccd_trial && i > cycles_/2 && ( numeric::random::rg().uniform() * cycles_ ) < i ) {
			ccd_trial->apply( pose );
		}
	}
}

void
LoopClosure::catch_fragment( Pose const& pose ) {
	bool success ( closure_frame_->steal( pose ) );
	runtime_assert( success );
}

void
LoopClosure::set_scorefxn( core::scoring::ScoreFunctionOP scorefxn ) {
	scorefxn_ = scorefxn;
	init_mc();
}

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols
