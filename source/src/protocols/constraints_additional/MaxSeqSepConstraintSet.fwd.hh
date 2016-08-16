// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/abinitio/MaxSeqSepConstraintSet.fwd.hh
/// @brief  forward declaration for MaxSeqSepConstraintSet class
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007


#ifndef INCLUDED_protocols_constraints_additional_MaxSeqSepConstraintSet_fwd_hh
#define INCLUDED_protocols_constraints_additional_MaxSeqSepConstraintSet_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace constraints_additional {

//forward declaration for private class
class MaxSeqSepConstraintSet;
typedef utility::pointer::shared_ptr< MaxSeqSepConstraintSet > MaxSeqSepConstraintSetOP;
typedef utility::pointer::shared_ptr< MaxSeqSepConstraintSet const> MaxSeqSepConstraintSetCOP;

}
}

#endif
