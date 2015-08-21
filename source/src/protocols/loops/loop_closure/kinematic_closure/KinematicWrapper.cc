// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/loop_closure/kinematic_closure/KinematicWrapper.cc
/// @brief KinematicWrapper methods implemented - this is a mover which simplifies use of KinematicMover loop modeling
/// @author Steven Lewis

// Unit Headers
#include <protocols/loops/loop_closure/kinematic_closure/KinematicWrapper.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>

// Package Headers

// Project Headers


#include <core/kinematics/MoveMap.hh>

// Utility Headers
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/exit.hh>
#include <numeric/random/random.hh>
//#include <utility/vector1.hh>

// option key includes
#include <basic/options/keys/loops.OptionKeys.gen.hh>

#include <protocols/loops/Loop.hh>
#include <utility/vector1.hh>

//Auto Headers


// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.loops.loop_closure.kinematic_closure.KinematicWrapper" );

namespace protocols {
namespace loops {
namespace loop_closure {
namespace kinematic_closure {

/// @details the responsiblity of apply() is to use the underlying
/// KinematicMover to close the loop, taking care of moving the KinematicMover's
/// begin/middle/end and ensuring that closure does occur.
void KinematicWrapper::apply( core::pose::Pose & pose ){

	core::Size counter(1);
	core::Size alc_start(0), alc_middle(0), alc_end(0); // three pivot residues for kinematic loop closure
	core::Size alc_end_in_vec(0), alc_start_in_vec(0); //alc_start's position in allowed_positions_ vector
	core::Size const npos(allowed_positions_.size());

	runtime_assert(!(npos < 3)); //need at least three res for pivots

	for ( /*counter*/; counter<=limit_; ++counter ) {
		//pick new points for kinematic closure
		// AS: the previous implementation had a "history bias" towards the N-terminus of the loop, as the start pivot can be anywhere between begin_loop and end_loop-2, while the choice of the end pivot depends on the start pivot
		// note that this implementation does not consider length restrictions for the loop, so it probably shouldn't be used for whole-protein ensemble generation -- this should be incorporated before putting the KinematicWrapper in charge of all pivot selections, including those from refine/LoopMover_KIC
		if ( basic::options::option[ basic::options::OptionKeys::loops::legacy_kic ]() || counter % 2 == 0 ) {
			alc_start_in_vec = numeric::random::rg().random_range(1,npos-2); //choose a random spot in allowed vector, but not near the end
			alc_start = allowed_positions_[alc_start_in_vec]; //store it as start
			alc_end_in_vec = numeric::random::rg().random_range(alc_start_in_vec+2, npos); //chose a post-start spot in allowed vector
			alc_end = allowed_positions_[alc_end_in_vec]; //store it as end
		} else {
			alc_end_in_vec = numeric::random::rg().random_range(3, npos); //chose a a random spot in allowed vector, but not near the start
			alc_end = allowed_positions_[alc_end_in_vec]; //store it as end
			alc_start_in_vec = numeric::random::rg().random_range(1,alc_end_in_vec-2); //choose a spot somewhere before the end position
			alc_start = allowed_positions_[alc_start_in_vec]; //store it as start
		}
		core::Size middle_offset = (alc_end - alc_start) / 2; //pick the natural middle between the two
		alc_middle = alc_start + middle_offset;


		//must guaruntee alc_middle is in vector
		bool ok(false);
		core::Size i(alc_start_in_vec+1);
		for ( ; i<alc_end_in_vec; ++i ) { //iterate through vector
			if ( alc_middle == allowed_positions_[i] ) { //if a match is found, no problem, we're done, set flag
				ok = true;
				break;
			} else if ( alc_middle < allowed_positions_[i] ) break; //if we overshot it, we have to fix it
		}
		if ( !ok ) { //if we didn't find an exact match, i points to the first element in allowed_positions_ larger than middle
			//there is an edge case which we have to handle first: if (i-1) == alc_start_in_vec (alc_start == alc_middle) then we have to choose alc_middle == allowed_positions[i]
			if ( alc_start_in_vec == (i-1) ) alc_middle = allowed_positions_[i];
			else {
				core::Size const diff1(allowed_positions_[i] - alc_middle), diff2(alc_middle - allowed_positions_[i-1]);
				alc_middle = (diff2 <= diff1 ? allowed_positions_[i-1] : allowed_positions_[i]);
				//take whatever vector element is closest to the calculated middle
			}
		}

		//alc_middle = allowed_positions_[RG.random_range(alc_start_in_vec+1, alc_end_in_vec-1)]; //?what would this do?

		//TR << "cycle " << counter << " pivots begin/middle/end "
		//  << alc_start << '/' << alc_middle << '/' << alc_end << std::endl;

		kinmover_->set_pivots(alc_start, alc_middle, alc_end);
		kinmover_->apply(pose);
		if ( kinmover_->last_move_succeeded() ) break; //we hope this will be the exit point
	}//end for many cycles

	if ( counter <= limit_ ) {
		TR << "KinematicMover closed in " << counter << " cycles; loop was begin/middle/end "
			<< alc_start << '/' << alc_middle << '/' << alc_end << std::endl;
		set_last_move_status( protocols::moves::MS_SUCCESS );
	} else if ( counter > limit_ ) {
		TR << "KinematicMover failed to close in " << limit_ << " cycles" << std::endl;
		set_last_move_status( protocols::moves::FAIL_RETRY );
	} else { utility_exit_with_message("How did we get here? - KinematicWrapper"); }

	return;
}//apply

std::string
KinematicWrapper::get_name() const {
	return "KinematicWrapper";
}

/// @details this function compares positions in the allowed_positions vector with the movemap, and removes non-mobile positions
void KinematicWrapper::respect_this_movemap( core::kinematics::MoveMapCOP mm )
{
	init_allowed_pos(); //reset the vector

	typedef utility::vector1<core::Size>::iterator iter;
	for ( iter i(allowed_positions_.begin()); i!= allowed_positions_.end(); ) {
		if ( mm->get_bb(*i) ) ++i;
		else i = allowed_positions_.erase(i);
	}

	TR << "respect_this_movemap has restricted loop pivots to these positions:";
	for ( iter i(allowed_positions_.begin()); i!= allowed_positions_.end(); ++i ) TR << " " << *i;
	TR << std::endl;

	return;
}

using namespace basic::options;
/// @brief ctor with Loop
KinematicWrapper::KinematicWrapper(
	KinematicMoverOP kinmover_in,
	protocols::loops::Loop loop_in,
	core::Size cycles
) : Mover(), kinmover_(kinmover_in), loop_begin_(loop_in.start()), loop_end_(loop_in.stop()),
	limit_( (cycles == 0) ? option[OptionKeys::loops::kinematic_wrapper_cycles].value() : cycles) //option or parameter
{
	ctor();
}

/// @brief ctor with explicit loop begin/end
KinematicWrapper::KinematicWrapper(
	KinematicMoverOP kinmover_in,
	core::Size loop_begin,
	core::Size loop_end,
	core::Size cycles
) : Mover(), kinmover_(kinmover_in), loop_begin_(loop_begin), loop_end_(loop_end),
	limit_( (cycles == 0) ? option[OptionKeys::loops::kinematic_wrapper_cycles].value() : cycles) //option or parameter
{
	ctor();
}

/// @details trivial wrapper around stuff needed in two ctors (de-duplicating the code)
void KinematicWrapper::ctor(){
	Mover::type( "KinematicWrapper" );
	runtime_assert((loop_end_ - loop_begin_ +1) >=3 );
	//should we clone (deep copy?) kinmover instead of just copying the OP?
	//NO, because that would complicate simulated annealing
	init_allowed_pos();
	return;
}

/// @details initializes the allowed_positions vector with every position in the loop as an allowed pivot
void KinematicWrapper::init_allowed_pos(){
	allowed_positions_.clear();
	allowed_positions_.reserve(loop_end_ - loop_begin_);
	for ( core::Size i(loop_begin_); i <= loop_end_; ++i ) allowed_positions_.push_back(i);
	return;
}

KinematicWrapper::~KinematicWrapper(){}

} // namespace kinematic_closure
} // namespace loop_closure
} // namespace loops
} // namespace protocols

