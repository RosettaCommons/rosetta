// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rigid/RigidBodyMover.fwd.hh
/// @brief  RigidBodyMover forward declarations header
/// @author

#ifndef INCLUDED_protocols_rigid_RigidBodyMover_fwd_hh
#define INCLUDED_protocols_rigid_RigidBodyMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace rigid {

//Forwards and OP typedefs
class RigidBodyMover;
typedef utility::pointer::shared_ptr< RigidBodyMover > RigidBodyMoverOP;
typedef utility::pointer::shared_ptr< RigidBodyMover const > RigidBodyMoverCOP;

class RigidBodyPerturbMover;
typedef utility::pointer::shared_ptr< RigidBodyPerturbMover > RigidBodyPerturbMoverOP;
typedef utility::pointer::shared_ptr< RigidBodyPerturbMover const > RigidBodyPerturbMoverCOP;

class RigidBodyPerturbRandomJumpMover;
typedef utility::pointer::shared_ptr< RigidBodyPerturbRandomJumpMover > RigidBodyPerturbRandomJumpMoverOP;
typedef utility::pointer::shared_ptr< RigidBodyPerturbRandomJumpMover const > RigidBodyPerturbRandomJumpMoverCOP;

class RigidBodyPerturbNoCenterMover;
typedef utility::pointer::shared_ptr< RigidBodyPerturbNoCenterMover > RigidBodyPerturbNoCenterMoverOP;
typedef utility::pointer::shared_ptr< RigidBodyPerturbNoCenterMover const > RigidBodyPerturbNoCenterMoverCOP;

class RigidBodyRandomizeMover;
typedef utility::pointer::shared_ptr< RigidBodyRandomizeMover > RigidBodyRandomizeMoverOP;
typedef utility::pointer::shared_ptr< RigidBodyRandomizeMover const > RigidBodyRandomizeMoverCOP;

class RigidBodySpinMover;
typedef utility::pointer::shared_ptr< RigidBodySpinMover > RigidBodySpinMoverOP;
typedef utility::pointer::shared_ptr< RigidBodySpinMover const > RigidBodySpinMoverCOP;

class RigidBodyDeterministicSpinMover;
typedef utility::pointer::shared_ptr< RigidBodySpinMover > RigidBodyDeterministicSpinMoverOP;
typedef utility::pointer::shared_ptr< RigidBodySpinMover const > RigidBodyDeterministicSpinMoverCOP;

class RigidBodyTiltMover;
typedef utility::pointer::shared_ptr< RigidBodyTiltMover > RigidBodyTiltMoverOP;
typedef utility::pointer::shared_ptr< RigidBodyTiltMover const > RigidBodyTiltMoverCOP;

class RigidBodyTransMover;
typedef utility::pointer::shared_ptr< RigidBodyTransMover > RigidBodyTransMoverOP;
typedef utility::pointer::shared_ptr< RigidBodyTransMover const > RigidBodyTransMoverCOP;

class UniformSphereTransMover;
typedef utility::pointer::shared_ptr< UniformSphereTransMover > UniformSphereTransMoverOP;
typedef utility::pointer::shared_ptr< UniformSphereTransMover const > UniformSphereTransMoverCOP;

class RigidBodyDofRandomizeMover;
typedef utility::pointer::shared_ptr< RigidBodyDofRandomizeMover > RigidBodyDofRandomizeMoverOP;
typedef utility::pointer::shared_ptr< RigidBodyDofRandomizeMover const > RigidBodyDofRandomizeMoverCOP;

class RigidBodyDofSeqRandomizeMover;
typedef utility::pointer::shared_ptr< RigidBodyDofSeqRandomizeMover > RigidBodyDofSeqRandomizeMoverOP;
typedef utility::pointer::shared_ptr< RigidBodyDofSeqRandomizeMover const > RigidBodyDofSeqRandomizeMoverCOP;

class RigidBodyDofTransMover;
typedef utility::pointer::shared_ptr< RigidBodyDofTransMover > RigidBodyDofTransMoverOP;
typedef utility::pointer::shared_ptr< RigidBodyDofTransMover const > RigidBodyDofTransMoverCOP;

class RigidBodyDofSeqTransMover;
typedef utility::pointer::shared_ptr< RigidBodyDofSeqTransMover > RigidBodyDofSeqTransMoverOP;
typedef utility::pointer::shared_ptr< RigidBodyDofSeqTransMover const > RigidBodyDofSeqTransMoverCOP;

class RigidBodyDofRandomTransMover;
typedef utility::pointer::shared_ptr< RigidBodyDofRandomTransMover > RigidBodyDofRandomTransMoverOP;
typedef utility::pointer::shared_ptr< RigidBodyDofRandomTransMover const > RigidBodyDofRandomTransMoverCOP;

class RigidBodyDofPerturbMover;
typedef utility::pointer::shared_ptr< RigidBodyDofPerturbMover > RigidBodyDofPerturbMoverOP;
typedef utility::pointer::shared_ptr< RigidBodyDofPerturbMover const > RigidBodyDofPerturbMoverCOP;

class RigidBodyDofSeqPerturbMover;
typedef utility::pointer::shared_ptr< RigidBodyDofSeqPerturbMover > RigidBodyDofSeqPerturbMoverOP;
typedef utility::pointer::shared_ptr< RigidBodyDofSeqPerturbMover const > RigidBodyDofSeqPerturbMoverCOP;

}//moves
}//protocols

#endif //INCLUDED_protocols_moves_RigidBodyMover_FWD_HH
