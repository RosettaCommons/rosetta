// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/constraints/SiteConstraintResidues.fwd.hh
/// @brief This class is an AmbiguousConstraint in which the set is comprised of AtomPairConstraints
/// @brief of an atom of interest in one chain versus the CA of of another sets of residues
/// @author Lei Shi (shilei@uw.edu)

#ifndef INCLUDED_core_scoring_constraints_SiteConstraintResidues_fwd_hh
#define INCLUDED_core_scoring_constraints_SiteConstraintResidues_fwd_hh

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace constraints {

class SiteConstraintResidues;

typedef utility::pointer::shared_ptr< SiteConstraintResidues > SiteConstraintResiduesOP;
typedef utility::pointer::shared_ptr< SiteConstraintResidues const > SiteConstraintResiduesCOP;

}
}
}

#endif
