// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief Cost computation for Gunn Moves
/// @author Oliver Lange

// Unit Headers
#include <protocols/simple_moves/SmoothFragmentMover.hh>

// Package Headers

// Project Headers
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/types.hh>
#include <basic/prof.hh>

#include <protocols/moves/Mover.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/random/random.hh>

//// C++ headers
#include <cstdlib>
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

/// @details Auto-generated virtual destructor
FragmentCost::~FragmentCost() {}


using namespace core;

SmoothFragmentMover::SmoothFragmentMover(
	core::fragment::FragSetCOP fragset,
	FragmentCostOP cost ) :
	ClassicFragmentMover( fragset, "SmoothFragmentMover_"+cost->type() ),
	cost_( cost )
{}


SmoothFragmentMover::SmoothFragmentMover(
	core::fragment::FragSetCOP fragset,
	core::kinematics::MoveMapCOP movemap,
	FragmentCostOP cost ) :
	ClassicFragmentMover( fragset, movemap, "SmoothFragmentMover_"+cost->type() ),
	cost_( cost )
{}

SmoothFragmentMover::SmoothFragmentMover(
	core::fragment::FragSetCOP fragset,
	core::kinematics::MoveMapCOP movemap,
	FragmentCostOP cost,
	std::string move_type ) :
	ClassicFragmentMover( fragset, movemap, move_type+"_"+cost->type() ),
	cost_( cost )
{}

std::string
SmoothFragmentMover::get_name() const {
	return "SmoothFragmentMover";
}

SmoothFragmentMover::~SmoothFragmentMover() {}

bool
SmoothFragmentMover::choose_fragment(
	core::fragment::FrameList const& frames,
	core::pose::Pose const& pose,
	Size &frame_num,
	Size &frag_num
) const
{

	PROF_START( basic::TEST4 );

	//std::cout << "SmoothFragmentMover::choose_fragment" << std::endl;
	typedef std::pair< Size, Size > FragID;
	utility::vector1< FragID > goodfrag;
	goodfrag.reserve( frames.size()*200 );

	Real costmin = 1000;
	FragID minfrag ( 0 , 0 );

	for ( Size fnr = 1; fnr <= frames.size(); fnr++ ) {
		//compute scores
		ScoreList scores;
		fragment::Frame const& frame( *( frames[ fnr ] ) );
		cost_->score( frame, pose, scores );

		for ( Size j = 1; j <= frame.nr_frags(); ++j ) {
			Real s = scores[ j ];
			//std::cout << "SmoothFragmentMover::choose_fragment: " << j << " score: " << s << std::endl;
			if ( s < costmin ) {
				costmin = s;
				minfrag = FragID( fnr, j );
			}
			if ( s < cost_->cutoff() ) {
				goodfrag.push_back( FragID( fnr, j ) );
			}
		}
		//choose randomly one fragment of all those that are below cutoff
		//or choose the best fragment. Fail if the minimal cost is > 12.0
	}

	PROF_STOP( basic::TEST4 );

	if ( goodfrag.size()< 1 ) {
		if ( costmin > 12. ) return false;
		frame_num = minfrag.first;
		frag_num = minfrag.second;
	} else {
		FragID choice = goodfrag[ static_cast< int >( numeric::random::rg().uniform() * goodfrag.size() )+1 ];
		frame_num = choice.first;
		frag_num = choice.second;
	}
	return true;
}


bool
SmoothFragmentMover::use_ss_length_screen() const
{
	return true;
}


} // simple_moves
} // protocols
