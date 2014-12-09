// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007
///


#ifndef INCLUDED_protocols_simple_moves_ConstraintFragmentMover_HH
#define INCLUDED_protocols_simple_moves_ConstraintFragmentMover_HH

// Unit Headers
//#include <protocols/simple_moves/ConstraintFragmentMover.fwd.hh>
#include <protocols/simple_moves/SmoothFragmentMover.hh>

// Package Headers
#include <protocols/simple_moves/FragmentMover.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/fragment/FragSet.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/MoveMap.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace simple_moves {

typedef utility::vector1< core::Real > ScoreList;


class ConstraintFragmentMover : public ClassicFragmentMover {
public:
	ConstraintFragmentMover(
		core::fragment::FragSetCOP fragset,
		core::kinematics::MoveMapCOP movemap,
		FragmentCostOP cost )
	:
		ClassicFragmentMover( fragset, movemap, "SmoothFragmentMover_"+cost->type() ),
		cost_( cost )
	{}

	//	 void apply( core::pose::Pose & );
protected:
	ConstraintFragmentMover(
					core::fragment::FragSetCOP fragset,
					core::kinematics::MoveMapCOP movemap,
					FragmentCostOP cost,
					std::string move_type ) :
		ClassicFragmentMover( fragset, movemap, move_type+"_"+cost->type() ),
		cost_( cost ) {}

	// frame and fragment of choice, returns false if no good fragment is found
	bool choose_fragment( core::fragment::FrameList const&, core::pose::Pose const&, Size &frame_num, Size &frag_num );

private:
	FragmentCostOP cost_;

	// choose randomly fragments that are below cutoff_
	core::Real cutoff_;

};

} //simple_moves
} //protocols

#endif
