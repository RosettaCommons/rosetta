// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file /rosetta/rosetta_source/src/devel/sewing/BridgeFragmentMover.hh
/// @brief Reads in fragments
/// @author Tim Jacobs

#ifndef BRIDGEFRAGMENTEVALUATOR_HH_
#define BRIDGEFRAGMENTEVALUATOR_HH_

#include <devel/sewing/BridgeFragmentMover.fwd.hh>

#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>
#include <core/fragment/FragSetCollection.hh>

namespace devel {
namespace sewing {

class BridgeFragmentMover : public protocols::moves::Mover{

public:

	BridgeFragmentMover(utility::vector1<core::fragment::FragSetOP> frag_sets);

	~BridgeFragmentMover();

	std::string get_name() const {
		return "BridgeFragmentMover";
	}

	void apply(core::pose::Pose &);

	core::pose::Pose constructClosedBundle(core::pose::Pose & pose, utility::vector1<core::pose::Pose> & loop_poses);

private:
	utility::vector1<core::fragment::FragSetOP> frag_sets_;

	//Numer of CA atoms we are using from each helix for RMSD calculations to bridge fragments
	core::Size num_helical_residues_;
	core::Size size_helical_window_;

};

} //sewing namespace
} //devel namespace

#endif /* BRIDGEFRAGMENTEVALUATOR_HH_ */
