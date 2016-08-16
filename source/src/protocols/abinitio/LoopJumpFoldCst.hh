// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file KinematicTaskCenter
/// @brief  this class will be handled to a SampleProtocol as a control instance
/// @details responsibilities:
///           know which chainbreaks to penalize and close
///           know which jumps to use during sampling, which (if any) to keep after loop-closing
///           supply a JumpMover if jumps should be moved
///           supply a MoveMap
///           supply a "StrictMoveMap": the protocol should not move anything that is dissallowed in strict_movemap(),
///                      it should try to move just stuff in movemap()
/// should this class also know how to ramp score terms ?
/// handle the titration of constraints ?
/// @author Oliver Lange


#ifndef INCLUDED_protocols_abinitio_LoopJumpFoldCst_hh
#define INCLUDED_protocols_abinitio_LoopJumpFoldCst_hh

// Unit Headers
#include <protocols/abinitio/LoopJumpFoldCst.fwd.hh>
#include <protocols/abinitio/KinematicTaskControl.hh>

// Package Headers
#include <protocols/jumping/JumpSetup.hh>
#include <protocols/abinitio/KinematicControl.hh>
#include <protocols/abinitio/Protocol.hh>
#include <protocols/abinitio/ResolutionSwitcher.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/loops/Loops.hh>
#include <core/fragment/SecondaryStructure.hh>

// ObjexxFCL Headers
//#include <ObjexxFCL/FArray1D.hh>
//#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <cstdlib>
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace abinitio {


class LoopJumpFoldCst : public KinematicTaskControl {
public:
	LoopJumpFoldCst(
		jumping::BaseJumpSetupOP jump_def,
		loops::Loops loops,
		ProtocolOP sampler,
		core::fragment::SecondaryStructureOP ss_def,
		core::Real coord_cst_weight,
		bool bCstAllAtom = false //false = CA atoms only
	) :
		KinematicTaskControl( sampler ),
		jump_def_ ( jump_def ),
		loops_( loops ),
		ss_def_( ss_def ),
		coordinate_constraint_weight_( coord_cst_weight ),
		bCstAllAtom_( bCstAllAtom )
	{
		dump_weights_file_ = "NO_DUMP";
	};

	~LoopJumpFoldCst();
	//@brief make a new KinematicControl...
	virtual KinematicControlOP new_kinematics( core::pose::Pose &pose );

	void set_coord_cst_weight_array( utility::vector1< core::Real > const& vec ) {
		coordinate_constraint_weights_ = vec;
	}

	void set_dump_weights_file( std::string const& str ) {
		dump_weights_file_ = str;
	}

	virtual std::string get_name() const;

protected:
	/// @brief heuristic to select subset of loops from loops_
	virtual void select_loops( loops::Loops& loops_select ) const;

	//@brief change fold-tree such that only loop regions move. keep jumps if they are not within the same rigid-region
	virtual bool add_rigidity_jumps( loops::Loops const& rigid, KinematicControlOP current_kinematics );


	virtual bool add_coord_cst(  loops::Loops const& loops, core::pose::Pose &pose );
	// virtual bool apply( core::pose::Pose &pose );

	virtual bool parse_jump_def( KinematicControlOP current_kinematics, core::kinematics::MoveMapOP );


private:
	jumping::BaseJumpSetupOP jump_def_;
	loops::Loops loops_; //if empty rebuild whole structure
	core::fragment::SecondaryStructureOP ss_def_;

protected:
	// a global weight
	core::Real coordinate_constraint_weight_;

	// a weight for each residue: if empty weight will be computed by heuristic
	utility::vector1< core::Real > coordinate_constraint_weights_;

private:
	// dump_weights to this file
	std::string dump_weights_file_;

	bool bCstAllAtom_;


};

}
}
#endif
