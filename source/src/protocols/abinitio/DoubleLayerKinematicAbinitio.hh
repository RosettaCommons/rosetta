// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file KinematicTaskCenter
/// @brief  this class will be handled to a SampleProtocol as a control instance
/// @details
/// this class will have two LOOP definitions:
/// loops ( defines which part has missing density and should start "extended"
/// rigids ( defines which part should be kept as a rigid (minimal) core )
/// the extended part will be sampled with stage1
/// following on this the non-rigid part will be sampled with stage2-stage4 and loop-closing

/// @author Oliver Lange


#ifndef INCLUDED_protocols_abinitio_DoubleLayerKinematicAbinitio_hh
#define INCLUDED_protocols_abinitio_DoubleLayerKinematicAbinitio_hh

// Unit Headers
#include <protocols/abinitio/DoubleLayerKinematicAbinitio.fwd.hh>
#include <protocols/abinitio/LoopJumpFoldCst.hh>

// Package Headers
#include <protocols/jumping/JumpSetup.hh>
#include <protocols/abinitio/KinematicControl.hh>
#include <protocols/abinitio/Protocol.hh>
#include <protocols/abinitio/ResolutionSwitcher.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/loops/Loops.hh>

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


class DoubleLayerKinematicAbinitio : public LoopJumpFoldCst {
public:
	DoubleLayerKinematicAbinitio(
		jumping::BaseJumpSetupOP jump_def,
		loops::Loops extended_loops,
		loops::Loops rigid_core,
		ProtocolOP sampler,
		ProtocolOP extended_chain_sampler,
		core::fragment::SecondaryStructureOP ss_def,
		core::Real coord_cst_weight,
		bool bCstAllAtom = false //false = CA atoms only
	) :
		LoopJumpFoldCst( jump_def, extended_loops, sampler, ss_def, coord_cst_weight, bCstAllAtom ),
		rigid_loops_( rigid_core ),
		extended_loops_( extended_loops ),
		stage1_sampler_( extended_chain_sampler )
	{}

	~DoubleLayerKinematicAbinitio();

	//@brief make a new KinematicControl...
	virtual KinematicControlOP new_kinematics( core::pose::Pose &pose );

	virtual std::string get_name() const;

	// virtual void init( core::pose::Pose const& pose );

protected:
	/// @brief heuristic to select subset of loops from loops_
	virtual void select_core_loops( loops::Loops& loops_select ) const;

	virtual bool inner_loop( core::pose::Pose& pose );
private:
	loops::Loops rigid_loops_; //if empty rebuild whole structure
	loops::Loops extended_loops_;
	ProtocolOP stage1_sampler_;
};

}
}
#endif
