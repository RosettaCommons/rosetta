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
///	  MultiResolutionProtocol
///   if your protocol wants to switch between full-atom and centroid representation derive it from this one
///   functionality to copy side-chains for fixed residues from an initial fa-pose to an  post-centroid fa-pose is provided
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_abinitio_ResolutionSwitcher_hh
#define INCLUDED_protocols_abinitio_ResolutionSwitcher_hh

// Unit Headers
#include <protocols/abinitio/ResolutionSwitcher.fwd.hh>

// Package Headers
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace abinitio {

class ResolutionSwitcher : public moves::Mover {
public:
	ResolutionSwitcher( core::pose::Pose const&, bool fullatom /*input pose*/, bool start_centroid, bool apply_to_centroid );
	~ResolutionSwitcher() {};

	//@brief if input was full-atom but we started (start_pose) from centroid we will copy side-chains
	// repacks all residues that have been moved between start_pose and pose
	virtual void apply( core::pose::Pose& );

	virtual std::string get_name() const;

	//@brief gives a starting pose ( with respect to setting in start_centroid )
	core::pose::Pose start_pose() const;

	core::pose::Pose const& init_pose() const {
		return init_pose_;
	}

	core::pose::Pose& init_pose() {
		return init_pose_;
	}

	void set_scorefxn( core::scoring::ScoreFunctionOP scorefxn_fa ) {
		scorefxn_fa_ = scorefxn_fa;
	}

	void set_map_cst_from_centroid_to_fa( bool const setting ){ map_cst_from_centroid_to_fa_ = setting; }


private:
	/// @brief true if apply() method is called on centroid pose
	bool apply_to_centroid_;


	/// @brief init_pose ( before the sampling started ) -- used to steal sidechains if fullatom
	core::pose::Pose init_pose_;

	/// @brief init_pose is full-atom
	bool init_fa_;

	/// @brief true if we want the start_pose() to be centroid
	bool start_centroid_;

	/// @brief full-atom scorefunction for repacking
	core::scoring::ScoreFunctionOP scorefxn_fa_;

	/// @brief repack sidechains at moved positions and x positions to left and right
	Size repack_buffer_;

	bool map_cst_from_centroid_to_fa_;

};

}
}

#endif

