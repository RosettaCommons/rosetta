// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/rigid/RigidBodyMover.fwd.hh
/// @brief  RigidBodyMover forward declarations header
/// @author

#ifndef INCLUDED_protocols_rigid_RigidBodyMover_fwd_hh
#define INCLUDED_protocols_rigid_RigidBodyMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols{
namespace rigid{

//Forwards and OP typedefs
class RigidBodyMover;
typedef utility::pointer::owning_ptr< RigidBodyMover > RigidBodyMoverOP;
typedef utility::pointer::owning_ptr< RigidBodyMover const > RigidBodyMoverCOP;

class RigidBodyPerturbMover;
typedef utility::pointer::owning_ptr< RigidBodyPerturbMover > RigidBodyPerturbMoverOP;
typedef utility::pointer::owning_ptr< RigidBodyPerturbMover const > RigidBodyPerturbMoverCOP;

class RigidBodyRandomTMHMover;
typedef utility::pointer::owning_ptr<RigidBodyRandomTMHMover> RigidBodyRandomTMHMoverOP;
typedef utility::pointer::owning_ptr<RigidBodyRandomTMHMover const> RigidBodyRandomTMHMoverCOP;

class RigidBodyPerturbRandomJumpMover;
typedef utility::pointer::owning_ptr< RigidBodyPerturbRandomJumpMover > RigidBodyPerturbRandomJumpMoverOP;
typedef utility::pointer::owning_ptr< RigidBodyPerturbRandomJumpMover const > RigidBodyPerturbRandomJumpMoverCOP;

class RigidBodyPerturbNoCenterMover;
typedef utility::pointer::owning_ptr< RigidBodyPerturbNoCenterMover > RigidBodyPerturbNoCenterMoverOP;
typedef utility::pointer::owning_ptr< RigidBodyPerturbNoCenterMover const > RigidBodyPerturbNoCenterMoverCOP;

class RigidBodyRandomizeMover;
typedef utility::pointer::owning_ptr< RigidBodyRandomizeMover > RigidBodyRandomizeMoverOP;
typedef utility::pointer::owning_ptr< RigidBodyRandomizeMover const > RigidBodyRandomizeMoverCOP;

class RigidBodySpinMover;
typedef utility::pointer::owning_ptr< RigidBodySpinMover > RigidBodySpinMoverOP;
typedef utility::pointer::owning_ptr< RigidBodySpinMover const > RigidBodySpinMoverCOP;

class RigidBodyTransMover;
typedef utility::pointer::owning_ptr< RigidBodyTransMover > RigidBodyTransMoverOP;
typedef utility::pointer::owning_ptr< RigidBodyTransMover const > RigidBodyTransMoverCOP;

class UniformSphereTransMover;
typedef utility::pointer::owning_ptr< UniformSphereTransMover > UniformSphereTransMoverOP;
typedef utility::pointer::owning_ptr< UniformSphereTransMover const > UniformSphereTransMoverCOP;

class RigidBodyDofRandomizeMover;
typedef utility::pointer::owning_ptr< RigidBodyDofRandomizeMover > RigidBodyDofRandomizeMoverOP;
typedef utility::pointer::owning_ptr< RigidBodyDofRandomizeMover const > RigidBodyDofRandomizeMoverCOP;

class RigidBodyDofSeqRandomizeMover;
typedef utility::pointer::owning_ptr< RigidBodyDofSeqRandomizeMover > RigidBodyDofSeqRandomizeMoverOP;
typedef utility::pointer::owning_ptr< RigidBodyDofSeqRandomizeMover const > RigidBodyDofSeqRandomizeMoverCOP;

class RigidBodyDofTransMover;
typedef utility::pointer::owning_ptr< RigidBodyDofTransMover > RigidBodyDofTransMoverOP;
typedef utility::pointer::owning_ptr< RigidBodyDofTransMover const > RigidBodyDofTransMoverCOP;

class RigidBodyDofSeqTransMover;
typedef utility::pointer::owning_ptr< RigidBodyDofSeqTransMover > RigidBodyDofSeqTransMoverOP;
typedef utility::pointer::owning_ptr< RigidBodyDofSeqTransMover const > RigidBodyDofSeqTransMoverCOP;

class RigidBodyDofRandomTransMover;
typedef utility::pointer::owning_ptr< RigidBodyDofRandomTransMover > RigidBodyDofRandomTransMoverOP;
typedef utility::pointer::owning_ptr< RigidBodyDofRandomTransMover const > RigidBodyDofRandomTransMoverCOP;

class RigidBodyDofPerturbMover;
typedef utility::pointer::owning_ptr< RigidBodyDofPerturbMover > RigidBodyDofPerturbMoverOP;
typedef utility::pointer::owning_ptr< RigidBodyDofPerturbMover const > RigidBodyDofPerturbMoverCOP;

class RigidBodyDofSeqPerturbMover;
typedef utility::pointer::owning_ptr< RigidBodyDofSeqPerturbMover > RigidBodyDofSeqPerturbMoverOP;
typedef utility::pointer::owning_ptr< RigidBodyDofSeqPerturbMover const > RigidBodyDofSeqPerturbMoverCOP;

}//moves
}//protocols

#endif //INCLUDED_protocols_moves_RigidBodyMover_FWD_HH
