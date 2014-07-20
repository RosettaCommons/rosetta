// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/MonteCarlo.tmpl.hh
/// @brief  implentation for MonteCarlo template functions
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_moves_MonteCarlo_tmpl_hh
#define INCLUDED_protocols_moves_MonteCarlo_tmpl_hh

// unit headers
#include <protocols/moves/MonteCarlo.hh>

// project headers
#include <core/pose/Pose.hh>

#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/PDB_Info.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MonteCarloStatus.hh>
#include <protocols/moves/MonteCarloExceptionConverge.fwd.hh>
#include <protocols/moves/MonteCarloExceptionConverge.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>



namespace protocols {
namespace moves {


/// @brief attach observer to last accepted conformation
/// @tparam ConformationObserver any class implementing <tt> void attach_to( Conformation & ) </tt>
template< typename ConformationObserver >
void
MonteCarlo::attach_observer_to_last_accepted_conformation( ConformationObserver & obs ) {
	obs.attach_to( last_accepted_pose_->conformation() );
}


/// @brief attach observer to lowest score conformation
/// @tparam ConformationObserver any class implementing <tt> void attach_to( Conformation & ) </tt>
template< typename ConformationObserver >
void
MonteCarlo::attach_observer_to_lowest_score_conformation( ConformationObserver & obs ) {
	obs.attach_to( lowest_score_pose_->conformation() );
}


/// @brief attach observer to last accepted pose
/// @tparam PoseObserver any class implementing <tt> void attach_to( Pose & ) </tt>
template< typename PoseObserver >
void
MonteCarlo::attach_observer_to_last_accepted_pose( PoseObserver & obs ) {
	obs.attach_to( *last_accepted_pose_ );
}


/// @brief attach observer to lowest score pose
/// @tparam PoseObserver any class implementing <tt> void attach_to( Pose & ) </tt>
template< typename PoseObserver >
void
MonteCarlo::attach_observer_to_lowest_score_pose( PoseObserver & obs ) {
	obs.attach_to( *lowest_score_pose_ );
}


} // namespace moves
} // protocols


#endif /* INCLUDED_protocols_moves_MonteCarlo_TMPL_HH */
