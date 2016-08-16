// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///	  Contains currently: LoopModeler
///
///
/// @author Vatsan Raman


#ifndef INCLUDED_devel_integrated_loop_CcdCloseMover_hh
#define INCLUDED_devel_integrated_loop_CcdCloseMover_hh

//Mover
#include <protocols/moves/Mover.hh>

//Core
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>

//LoopModel & Loops
#include <devel/integrated_loop/LoopModeler.hh>
#include <devel/integrated_loop/LoopManager.hh>

// protocols header
#include <protocols/loops/looprelax_protocols.hh>


// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/vector1.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <vector>

using core::Real;

namespace protocols {
namespace moves {

class CcdCloseMover;
typedef utility::pointer::owning_ptr< protocols::moves::CcdCloseMover >  CcdCloseMoverOP;

class CcdCloseMover: public protocols::moves::Mover {
public:
	CcdCloseMover(){}

 CcdCloseMover(
	 core::pose::Pose pose,
	 core::kinematics::MoveMap movemap,
	 protocols::Loop ThisLoop
 ):Mover(),
	 pose_( pose ),
	 movemap_( movemap ),
	 ThisLoop_( ThisLoop )
	{
		Mover::type("CcdCloseMover");
		set_default_param();
	}


	void set_default_param()
	{
	// param for ccd_closure
	ccd_cycles = 100;
	ccd_tol = 0.01;
	rama_check = true;
	max_rama_score_increase = 2.0;
	max_total_delta_helix = 10.0;
	max_total_delta_strand = 50.0;
	max_total_delta_loop = 75.0;
	}


	void apply(
		core::pose::Pose & pose_
	);

	private:
	core::pose::Pose pose_;
	core::kinematics::MoveMap movemap_;
	protocols::Loop ThisLoop_;

	// param for ccd_closure
	Size ccd_cycles; // num of cycles of ccd_moves
	Real ccd_tol; // criterion for a closed loop
	bool rama_check;
	Real max_rama_score_increase; // dummy number when rama_check is false
	Real max_total_delta_helix; // max overall angle changes for a helical residue
	Real max_total_delta_strand; // ... for a residue in strand
	Real max_total_delta_loop; // ... for a residue in loop


	// output for ccd_closure
	//	Real forward_deviation, backward_deviation; // actually loop closure msd, both dirs
	//	Real torsion_delta, rama_delta; // actually torsion and rama score changes, averaged by loop_size


};

} //moves
} //protocols

#endif
