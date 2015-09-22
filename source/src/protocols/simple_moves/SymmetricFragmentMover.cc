// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Inserts a Fragment into a Pose, similar to old Rosetta++ main_frag_trial algorithm.
/// @author Oliver Lange

// Unit Headers
#include <protocols/simple_moves/SymmetricFragmentMover.hh>

// Package Headers

// Project Headers
#include <core/fragment/Frame.hh>

// Utility headers
#include <numeric/random/random.hh>
#include <basic/Tracer.hh>

#include <core/fragment/FragData.hh>
#include <utility/vector1.hh>


// ObjexxFCL Headers

// C++ headers

namespace protocols {
namespace simple_moves {


using namespace core;
using namespace fragment;
using namespace basic;

static THREAD_LOCAL basic::Tracer tr( "protocols.simple_moves.FragmentMover" );

std::string
SymmetricFragmentMover::get_name() const {
	return "SymmetricFragmentMover";
}

bool
SymmetricFragmentMover::apply_fragment(
	core::fragment::Frame const& frame,
	Size frag_num,
	core::kinematics::MoveMap const& movemap,
	core::pose::Pose &pose
) const {
	bool success = ClassicFragmentMover::apply_fragment( frame, frag_num, movemap, pose );
	if ( success ) {
		Size new_start( 0 );
		if ( frame.start() >= image_start_ ) {
			new_start = frame.start() - image_start_ + 1;
		} else {
			new_start = frame.start() + image_start_ - 1;
		}
		if ( frame.is_continuous() ) {
			frame.fragment( frag_num ).apply( movemap, pose, new_start, new_start + frame.length() - 1 );
		} else {
			tr.Warning << "WARNING: symmetric mover did not copy fragment move for non-continous fragment, only applied on monomer" << std::endl;
		}
	}
	return success;
}


} // simple_moves
} // protocols
