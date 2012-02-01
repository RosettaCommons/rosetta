// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   
/// @brief  
/// @author Javier Castellanos ( javiercv@uw.edu )


// Package Header
#include <devel/constrained_sequence_design/SequenceConstraintSet.hh>

// Project Header
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

// STL Header
#include <set>

namespace devel {
namespace constrained_sequence_design {

SequenceConstraintSet::SequenceConstraintSet()
{ }

SequenceConstraintSet::~SequenceConstraintSet()
{ }

SequenceConstraintSet::SequenceConstraintSet(const SequenceConstraintSet& cs) 
	: constraints_( cs.constraints_ ) 
{ }

typedef utility::pointer::owning_ptr< SequenceConstraintSet > SequenceConstraintSetOP;
typedef utility::pointer::owning_ptr< SequenceConstraintSet const > SequenceConstraintSetCOP;

} // namespace devel
} // namespace constrained_sequence_design
