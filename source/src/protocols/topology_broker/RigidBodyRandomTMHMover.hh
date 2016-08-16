// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Stephanie and Sam DeLuca


#ifndef INCLUDED_protocols_topology_broker_RigidBodyRandomTMHMover_hh
#define INCLUDED_protocols_topology_broker_RigidBodyRandomTMHMover_hh

// Unit headers
#include <protocols/topology_broker/RigidBodyRandomTMHMover.fwd.hh>
#include <protocols/topology_broker/TopologyClaimer.fwd.hh>
#include <protocols/rigid/RigidBodyMover.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/conformation/symmetry/SymDof.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <map>

namespace protocols {
namespace topology_broker {

class RigidBodyRandomTMHMover : public rigid::RigidBodyMover{
public:
	RigidBodyRandomTMHMover();
	RigidBodyRandomTMHMover(core::Real max_trans, core::Real rotation_mag, core::Real translation_mag, core::Size tmhelix,
		protocols::topology_broker::TopologyClaimerOP claimer);
	void apply(core::pose::Pose& pose);
	virtual std::string get_name() const;
	virtual ~RigidBodyRandomTMHMover();

private:

	core::Real max_trans_;
	core::Real trans_mag_in_;
	core::Real rot_mag_in_;
	core::Size num_jump_;
	topology_broker::TopologyClaimerOP claimer_;

};

}
}

#endif
