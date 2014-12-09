// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /rosetta/rosetta_source/src/devel/sewing/BridgeFragmentMover.fwd.hh
/// @brief 
/// @author Tim Jacobs

#ifndef BRIDGEFRAGMENTMOVER_FWD_HH_
#define BRIDGEFRAGMENTMOVER_FWD_HH_

#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace sewing {

class BridgeFragmentMover;
typedef utility::pointer::owning_ptr< BridgeFragmentMover > BridgeFragmentMoverOP;
typedef utility::pointer::owning_ptr< BridgeFragmentMover const > BridgeFragmentMoverCOP;

} //sewing namespace
} //devel namespace

#endif /* BRIDGEFRAGMENTMOVER_FWD_HH_ */
