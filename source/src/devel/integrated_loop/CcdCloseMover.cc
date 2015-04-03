// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///	  Contains currently: LoopModeler
///
///
/// @author Vatsan Raman
#include <devel/integrated_loop/CcdCloseMover.hh>
#include <core/pose/Pose.hh>
#include <protocols/loops/ccd_closure.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>

using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;

namespace protocols {
namespace moves {
void CcdCloseMover::apply(
	core::pose::Pose & pose_
)

{


	// output for ccd_closure
	Real forward_deviation, backward_deviation; // actually loop closure msd, both dirs
	Real torsion_delta, rama_delta; // actually torsion and rama score changes, averaged by loop

	Size loop_begin( ThisLoop_.loop_begin() );
	Size loop_end( ThisLoop_.loop_end() );
	Size cutpoint( ThisLoop_.cutpoint() );

	// ccd close this loop
	protocols::loops::fast_ccd_loop_closure( pose_, movemap_, loop_begin, loop_end, cutpoint,
		ccd_cycles, ccd_tol, rama_check, max_rama_score_increase, max_total_delta_helix,
			max_total_delta_strand, max_total_delta_loop, forward_deviation,
			backward_deviation, torsion_delta, rama_delta );
}


}//moves
}//protocols
