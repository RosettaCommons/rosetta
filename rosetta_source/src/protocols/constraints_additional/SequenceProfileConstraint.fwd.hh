// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/constraints_additional/SequenceProfileConstraint.fwd.hh
/// @author ashworth


#ifndef INCLUDED_protocols_constraints_additional_SequenceProfileConstraint_fwd_hh
#define INCLUDED_protocols_constraints_additional_SequenceProfileConstraint_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace constraints_additional {

class SequenceProfileConstraint; // fwd declaration
typedef utility::pointer::owning_ptr< SequenceProfileConstraint > SequenceProfileConstraintOP;
typedef utility::pointer::owning_ptr< SequenceProfileConstraint const > SequenceProfileConstraintCOP;

} // namespace constraints_additional
} // namespace protocols

#endif // INCLUDED_protocols_constraints_additional_SequenceProfileConstraint_FWD_HH
