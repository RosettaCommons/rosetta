// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @author Oliver Lange
///
/// this class replaces the frag_close routines in jumping_pairings.cc
/// the short loop is copied into a special purpose pose that just contains the loop-fragment
/// with the hope that things are speeded up... which may or may not be true!
/// only the linear chainbreak score is used.

#ifndef INCLUDED_protocols_loops_loop_closure_ccd_ShortLoopClosure_hh
#define INCLUDED_protocols_loops_loop_closure_ccd_ShortLoopClosure_hh

// Unit Headers
#include <protocols/loops/loop_closure/ccd/ShortLoopClosure.fwd.hh>

// Package Headers
#include <protocols/loops/loop_closure/ccd/LoopClosure.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

class ShortLoopClosure : public LoopClosure {
	typedef LoopClosure Parent;
public:

	//c'stor
	ShortLoopClosure(
		core::fragment::FragSetCOP fragset,
		Loop loop_def,
		core::kinematics::MoveMapCOP movemap
	);

	//@brief run protocol on pose
	virtual bool apply( core::pose::Pose const& pose );

	/// @brief save the loop-fragment in closure_frames_
	// overwritten to realign frames to target sequence
	virtual void catch_fragment( core::pose::Pose const& short_pose );
private:
	Loop orig_loop_;
};

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols

#endif //INCLUDED_protocols_loops_loop_closure_ccd_ShortLoopClosure_hh

