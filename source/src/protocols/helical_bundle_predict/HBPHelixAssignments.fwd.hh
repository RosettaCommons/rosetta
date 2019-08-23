// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HBPHelixAssignments.fwd.hh
/// @brief A class for storing the helix assignments for a pose.  This can represent those proposed from an input file, or those at the current state of a trajectory.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_helical_bundle_predict_HBPHelixAssignments_fwd_hh
#define INCLUDED_protocols_helical_bundle_predict_HBPHelixAssignments_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

class HBPHelixAssignmentsTests; //Forward declaration needed for friendship.

// Forward
namespace protocols {
namespace helical_bundle_predict {

class HBPHelixParameters;
class HBPHelix;
class HBPHelixAssignments;

typedef utility::pointer::shared_ptr< HBPHelixParameters > HBPHelixParametersOP;
typedef utility::pointer::shared_ptr< HBPHelixParameters const > HBPHelixParametersCOP;

typedef utility::pointer::shared_ptr< HBPHelix > HBPHelixOP;
typedef utility::pointer::shared_ptr< HBPHelix const > HBPHelixCOP;

typedef utility::pointer::shared_ptr< HBPHelixAssignments > HBPHelixAssignmentsOP;
typedef utility::pointer::shared_ptr< HBPHelixAssignments const > HBPHelixAssignmentsCOP;

} //protocols
} //helical_bundle_predict

#endif //INCLUDED_protocols_helical_bundle_predict_HBPHelixAssignments_fwd_hh
