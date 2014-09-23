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
/// @detailed responsibilities:
///           know which chainbreaks to penalize and close
///           know which jumps to use during sampling, which (if any) to keep after loop-closing
///           supply a JumpMover if jumps should be moved
///           supply a MoveMap
///           supply a "StrictMoveMap": the protocol should not move anything that is dissallowed in strict_movemap(),
///                      it should try to move just stuff in movemap()
/// should this class also know how to ramp score terms ?
/// handle the titration of constraints ?
/// @author Oliver Lange


#ifndef INCLUDED_protocols_abinitio_KinematicTaskControl_hh
#define INCLUDED_protocols_abinitio_KinematicTaskControl_hh

// Unit Headers
#include <protocols/abinitio/KinematicTaskControl.fwd.hh>

// Package Headers
// AUTO-REMOVED #include <protocols/jumping/JumpSetup.hh>
#include <protocols/abinitio/KinematicControl.hh>
#include <protocols/abinitio/Protocol.hh>
#include <protocols/checkpoint/CheckPointer.hh>
#include <protocols/abinitio/ResolutionSwitcher.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// AUTO-REMOVED #include <protocols/loops/Loops.hh>

// ObjexxFCL Headers
//#include <ObjexxFCL/FArray1D.hh>
//#include <ObjexxFCL/FArray2D.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
// AUTO-REMOVED #include <cstdlib>
#include <string>

#include <protocols/loops/Loops.fwd.hh>
#include <utility/vector1.hh>



namespace protocols {
namespace abinitio {

class KinematicTaskControl : public Protocol {
public:
	KinematicTaskControl() :
		current_kinematics_( NULL ),
		sampling_protocol_( NULL )
	{
		return_centroid( false ); //default
	}

	KinematicTaskControl( ProtocolOP sampler ) :
		sampling_protocol_( sampler )
	{
		return_centroid( false );//default
	}

	~KinematicTaskControl();

	//@brief generate a new KinematicControl object
	virtual KinematicControlOP new_kinematics( core::pose::Pose &pose ) = 0;

	virtual void apply( core::pose::Pose &pose );
	virtual std::string get_name() const;

	virtual void init( core::pose::Pose const& pose );

  KinematicControlOP current_kinematics() {
    runtime_assert( current_kinematics_ != 0 );
    return current_kinematics_;
  }

	void set_input_pose_is_fa( bool setting = true ) {
		b_input_is_fullatom_ = setting;
	}

// 	virtual StructureStore const& structure_store() const {
// 		return sampling_protocol_->structure_store();
// 	}

// 	virtual StructureStore& structure_store() {
// 		return sampling_protocol_->structure_store();
// 	}

	ProtocolOP sampling_protocol() {
		return sampling_protocol_;
	}

	virtual checkpoint::CheckPointer &get_checkpoints() { return sampling_protocol_->get_checkpoints(); }

protected:
	virtual bool inner_loop( core::pose::Pose &pose );
	mutable KinematicControlOP current_kinematics_;
	void 	set_extended_torsions_and_idealize_loops( core::pose::Pose& pose, loops::Loops loops ) const;
	ResolutionSwitcher& res_switch() {
		return *res_switch_;
	}

private:
	ResolutionSwitcherOP res_switch_;
	ProtocolOP sampling_protocol_;
	bool b_input_is_fullatom_;
};



}
}

#endif
