// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_rbsegment_relax_RBSegment_fwd_hh
#define INCLUDED_protocols_rbsegment_relax_RBSegment_fwd_hh


// Rosetta Headers

// ObjexxFCL Headers

// Utility Headers
#include <utility/pointer/owning_ptr.fwd.hh>

#include <utility/vector1_bool.hh>


// C++ Headers

namespace protocols {
namespace rbsegment_relax {

class RBSegmentMover;
typedef utility::pointer::shared_ptr< RBSegmentMover >  RBSegmentMoverOP;

class RBSegment;
typedef utility::pointer::shared_ptr< RBSegment >  RBSegmentOP;
typedef utility::vector1< RBSegment > CompositeSegment;
typedef utility::vector1< RBSegment >::iterator RBIt;
typedef utility::vector1< RBSegment >::const_iterator RBConsIt;

class RBSegmentRelax;


} // ns moves
} // ns protocols

#endif
