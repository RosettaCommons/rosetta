// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rigid_body/RigidBodyUtil.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rotamer_sampler_rigid_body_RigidBodyUtil_HH
#define INCLUDED_protocols_rotamer_sampler_rigid_body_RigidBodyUtil_HH

#include <core/id/AtomID.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/types.hh>
#include <utility>

using namespace core;

namespace protocols {
namespace rotamer_sampler {
namespace rigid_body {

	void
	get_atom_coordinates( utility::vector1< std::pair < id::AtomID, numeric::xyzVector< core::Real > > > & xyz_list,
												Size const & seq_num,
												core::conformation::Residue const & rsd_at_origin,
												core::kinematics::Stub const & moving_res_base_stub );

} //rigid_body
} //rotamer_sampler
} //protocols

#endif
