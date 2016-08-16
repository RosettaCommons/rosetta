// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DisulfJumpClaimer
/// @brief  Claimer for disulfide jump sampling
/// @details responsibilities:
/// @author Robert Vernon


#ifndef INCLUDED_protocols_topology_broker_DisulfJumpClaimer_hh
#define INCLUDED_protocols_topology_broker_DisulfJumpClaimer_hh


// Unit Headers
#include <protocols/topology_broker/DisulfJumpClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/FragmentClaimer.hh>
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>
#include <protocols/topology_broker/claims/JumpClaim.fwd.hh>
#include <protocols/topology_broker/weights/AbinitioMoverWeight.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/jumping/DisulfPairingsList.fwd.hh>
#include <protocols/jumping/DisulfPairingLibrary.hh>

#ifdef WIN32
#include <protocols/jumping/DisulfPairingsList.hh>
#endif

#include <core/fragment/FrameList.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/vector1.fwd.hh>

//C++ Headers
#include <string>

namespace protocols {
namespace topology_broker {

/// @brief Claimer that works with the old system of BaseJumpSetup
/// it supports only JumpFrames of type  [ BBTorsion ] UpJump DownJump [ BBTorsion ]
/// the class JumpSample is still used to transport the information jumps and jump_atoms, but cuts are ignored
/// all functionality of JumpSample is not used anymore
class DisulfJumpClaimer : public FragmentClaimer {
	typedef FragmentClaimer Parent;
public:
	DisulfJumpClaimer(); //for factory
	~DisulfJumpClaimer();

	virtual TopologyClaimerOP clone() const;

	void generate_jump_frags( protocols::jumping::DisulfPairingLibrary const& lib, core::fragment::FrameList& all_frames) const;

	virtual void generate_claims( claims::DofClaims& );

	virtual void new_decoy( core::pose::Pose const& );
	virtual void new_decoy();

	/// @brief type() is specifying the output name of the TopologyClaimer
	virtual std::string type() const {
		return _static_type_name();
	}

	virtual void initialize_dofs( core::pose::Pose&, claims::DofClaims const& init_claims, claims::DofClaims& failed_to_init );

	static std::string _static_type_name() {
		return "DisulfJumpClaimer";
	}

	virtual bool read_tag( std::string tag, std::istream& is );

private:
	utility::vector1< claims::JumpClaimOP > local_disulf_data_;
	utility::vector1< protocols::jumping::DisulfPairing > all_jump_pairings_;
	core::fragment::FrameList all_frames_;

	bool bKeepJumpsFromInputPose_;
}; //class DisulfJumpClaimer

}
}

#endif
