// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


#ifndef INCLUDED_core_scoring_constraints_DOF_Constraint_fwd_hh
#define INCLUDED_core_scoring_constraints_DOF_Constraint_fwd_hh

// Utility header
#include <utility/pointer/owning_ptr.hh>

#include <utility/vector1.fwd.hh>


namespace core {
namespace scoring {
namespace constraints {

class DOF_Constraint;

typedef utility::pointer::shared_ptr< DOF_Constraint > DOF_ConstraintOP;
typedef utility::pointer::shared_ptr< DOF_Constraint const > DOF_ConstraintCOP;
typedef utility::vector1< DOF_ConstraintOP > DOF_ConstraintOPs;
typedef utility::vector1< DOF_ConstraintCOP > DOF_ConstraintCOPs;

}
}
}

#endif
