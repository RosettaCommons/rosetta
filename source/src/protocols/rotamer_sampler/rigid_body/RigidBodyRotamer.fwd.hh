// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rigid_body/RigidBodyRotamer.fwd.hh
/// @brief 
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rotamer_sampler_rigid_body_RigidBodyRotamer_FWD_HH
#define INCLUDED_protocols_rotamer_sampler_rigid_body_RigidBodyRotamer_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace rotamer_sampler {
namespace rigid_body {
	
	class RigidBodyRotamer;
	typedef utility::pointer::owning_ptr< RigidBodyRotamer > RigidBodyRotamerOP;
	typedef utility::pointer::owning_ptr< RigidBodyRotamer const > RigidBodyRotamerCOP;
	
} //rigid_body 
} //rotamer_sampler 
} //protocols 

#endif
