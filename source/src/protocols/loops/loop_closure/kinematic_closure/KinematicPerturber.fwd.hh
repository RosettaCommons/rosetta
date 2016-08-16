// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.fwd.hh
/// @brief  KinematicPerturber forward declarations header
/// @author Florian Richter (floric@u.washington.edu), march 2009
/// @author Amelie Stein (amelie.stein@ucsf.edu), Oct 2012


#ifndef INCLUDED_protocols_loops_loop_closure_kinematic_closure_KinematicPerturber_fwd_hh
#define INCLUDED_protocols_loops_loop_closure_kinematic_closure_KinematicPerturber_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace loops {
namespace loop_closure {
namespace kinematic_closure {

//Forwards and OP typedefs
class KinematicPerturber;
typedef utility::pointer::shared_ptr< KinematicPerturber > KinematicPerturberOP;
typedef utility::pointer::shared_ptr< KinematicPerturber const > KinematicPerturberCOP;

class TorsionSamplingKinematicPerturber;
typedef utility::pointer::shared_ptr< TorsionSamplingKinematicPerturber > TorsionSamplingKinematicPerturberOP;
typedef utility::pointer::shared_ptr< TorsionSamplingKinematicPerturber const > TorsionSamplingKinematicPerturberCOP;

class SmallMoveKinematicPerturber;
typedef utility::pointer::shared_ptr< SmallMoveKinematicPerturber > SmallMoveKinematicPerturberOP;
typedef utility::pointer::shared_ptr< SmallMoveKinematicPerturber const > SmallMoveKinematicPerturberCOP;

class PivotBasedKinematicPerturber;
typedef utility::pointer::shared_ptr< PivotBasedKinematicPerturber > PivotBasedKinematicPerturberOP;
typedef utility::pointer::shared_ptr< PivotBasedKinematicPerturber const > PivotBasedKinematicPerturberCOP;

class VicinitySamplingKinematicPerturber;
typedef utility::pointer::shared_ptr< VicinitySamplingKinematicPerturber > VicinitySamplingKinematicPerturberOP;
typedef utility::pointer::shared_ptr< VicinitySamplingKinematicPerturber const > VicinitySamplingKinematicPerturberCOP;

class TorsionSweepingKinematicPerturber;
typedef utility::pointer::shared_ptr< TorsionSweepingKinematicPerturber > TorsionSweepingKinematicPerturberOP;
typedef utility::pointer::shared_ptr< TorsionSweepingKinematicPerturber const > TorsionSweepingKinematicPerturberCOP;

class NeighborDependentTorsionSamplingKinematicPerturber;
typedef utility::pointer::shared_ptr< NeighborDependentTorsionSamplingKinematicPerturber > NeighborDependentTorsionSamplingKinematicPerturberOP;
typedef utility::pointer::shared_ptr< NeighborDependentTorsionSamplingKinematicPerturber const > NeighborDependentTorsionSamplingKinematicPerturberCOP;

class NeighborDependentTabooSamplingKinematicPerturber;
typedef utility::pointer::shared_ptr< NeighborDependentTabooSamplingKinematicPerturber > NeighborDependentTabooSamplingKinematicPerturberOP;
typedef utility::pointer::shared_ptr< NeighborDependentTabooSamplingKinematicPerturber const > NeighborDependentTabooSamplingKinematicPerturberCOP;

class TorsionRestrictedKinematicPerturber;
typedef utility::pointer::shared_ptr< TorsionRestrictedKinematicPerturber > TorsionRestrictedKinematicPerturberOP;
typedef utility::pointer::shared_ptr< TorsionRestrictedKinematicPerturber const > TorsionRestrictedKinematicPerturberCOP;

class TabooSamplingKinematicPerturber;
typedef utility::pointer::shared_ptr< TabooSamplingKinematicPerturber > TabooSamplingKinematicPerturberOP;
typedef utility::pointer::shared_ptr< TabooSamplingKinematicPerturber const > TabooSamplingKinematicPerturberCOP;

} // namespace kinematic_closure
} // namespace loop_closure
} // namespace loops
} // namespace protocols

#endif //INCLUDED_protocols_loops_loop_closure_kinematic_closure_KinematicPerturber_FWD_HH
