// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/kinematic_closure/KinematicPerturber.fwd.hh
/// @brief  KinematicPerturber forward declarations header
/// @author Florian Richter (floric@u.washington.edu), march 2009


#ifndef INCLUDED_protocols_moves_kinematic_closure_KinematicPerturber_fwd_hh
#define INCLUDED_protocols_moves_kinematic_closure_KinematicPerturber_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols{
namespace moves{
namespace kinematic_closure{

//Forwards and OP typedefs
class KinematicPerturber;
typedef utility::pointer::owning_ptr< KinematicPerturber > KinematicPerturberOP;
typedef utility::pointer::owning_ptr< KinematicPerturber const > KinematicPerturberCOP;

class TorsionSamplingKinematicPerturber;
typedef utility::pointer::owning_ptr< TorsionSamplingKinematicPerturber > TorsionSamplingKinematicPerturberOP;
typedef utility::pointer::owning_ptr< TorsionSamplingKinematicPerturber const > TorsionSamplingKinematicPerturberCOP;

class TorsionSweepingKinematicPerturber;
typedef utility::pointer::owning_ptr< TorsionSweepingKinematicPerturber > TorsionSweepingKinematicPerturberOP;
typedef utility::pointer::owning_ptr< TorsionSweepingKinematicPerturber const > TorsionSweepingKinematicPerturberCOP;


}//kinematic_closure
}//moves
}//protocols

#endif //INCLUDED_protocols_moves_kinematic_closure_KinematicPerturber_FWD_HH
