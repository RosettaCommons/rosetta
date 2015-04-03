// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @details responsibilities:
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_JumpClaimer_hh
#define INCLUDED_protocols_topology_broker_JumpClaimer_hh


// Unit Headers
#include <protocols/topology_broker/FragmentJumpClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/FragmentClaimer.hh>
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>
#include <protocols/topology_broker/weights/AbinitioMoverWeight.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <protocols/jumping/JumpSetup.hh>
#include <protocols/jumping/JumpSample.hh>

// ObjexxFCL Headers

// Utility headers
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
//#include <basic/Tracer.hh>
//#include <basic/options/option.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


//#include <basic/options/option_macros.hh>

//// C++ headers
//#include <fstream>


// option key includes


namespace protocols {
namespace topology_broker {

/// @brief Claimer that works with the old system of BaseJumpSetup
/// it supports only JumpFrames of type  [ BBTorsion ] UpJump DownJump [ BBTorsion ]
/// the class JumpSample is still used to transport the information jumps and jump_atoms, but cuts are ignored
/// all functionality of JumpSample is not used anymore
class FragmentJumpClaimer : public FragmentClaimer {
	typedef FragmentClaimer Parent;
public:
	FragmentJumpClaimer(); //for factory
	FragmentJumpClaimer( jumping::BaseJumpSetupOP jump_def, std::string const& mover_tag = "JumpMove",  weights::AbinitioMoverWeightOP weight = NULL );
	FragmentJumpClaimer( FragmentJumpClaimer const & src );
	~FragmentJumpClaimer();

	virtual TopologyClaimerOP clone() const;

	virtual void generate_claims( claims::DofClaims& );

	virtual void new_decoy( core::pose::Pose const& );
	virtual void new_decoy();

	/// @brief type() is specifying the output name of the TopologyClaimer
	virtual std::string type() const {
		return _static_type_name();
	}

	virtual void initialize_dofs( core::pose::Pose&, claims::DofClaims const& init_claims, claims::DofClaims& failed_to_init );

	static std::string _static_type_name() {
		return "FragmentJumpClaimer";
	}

	virtual bool read_tag( std::string tag, std::istream& is );

protected:
	void set_jump_def( jumping::BaseJumpSetupOP jump_def ) {
		jump_def_ = jump_def;
	}

	virtual void generate_claims( claims::DofClaims&, std::string, std::string );

	jumping::BaseJumpSetupOP jump_def() {
		return jump_def_;
	}

	void set_keep_jumps_from_input_pose( bool setting ) {
		bKeepJumpsFromInputPose_ = setting;
	}

	void init_jumps();

private:
	jumping::BaseJumpSetupOP jump_def_;
	jumping::JumpSample current_jumps_;
	simple_moves::ClassicFragmentMoverOP init_mover_;
	bool bKeepJumpsFromInputPose_;
	core::pose::Pose input_pose_;
	bool discard_jumps_;
}; //class JumpClaimer


}
}

#endif
