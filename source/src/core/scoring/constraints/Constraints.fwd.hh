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


#ifndef INCLUDED_core_scoring_constraints_Constraints_fwd_hh
#define INCLUDED_core_scoring_constraints_Constraints_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace constraints {

class Constraints;
typedef utility::pointer::shared_ptr< Constraints > ConstraintsOP;
typedef utility::pointer::shared_ptr< Constraints const > ConstraintsCOP;


} // constraints
} // scoring
} // core

#endif
