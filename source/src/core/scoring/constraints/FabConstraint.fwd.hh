// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/constraints/FabConstraint.fwd.hh
/// @brief This class is specific to antibodies and penalizes presence of residues flanking
/// @brief antibody cdr residues at Antigen-Antibody interfaces (ported from Fab constraint
/// @brief in rosetta++ which uses a constant constraint score of 0.5/flanking residue)
/// @author Krishna Kilambi (kkpraneeth@jhu.edu, April 2012)

#ifndef INCLUDED_core_scoring_constraints_FabConstraint_fwd_hh
#define INCLUDED_core_scoring_constraints_FabConstraint_fwd_hh

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace core {
    namespace scoring {
        namespace constraints {
            
            class FabConstraint;
            
            typedef utility::pointer::shared_ptr< FabConstraint > FabConstraintOP;
            typedef utility::pointer::shared_ptr< FabConstraint const > FabConstraintCOP;
            
        }
    }
}

#endif
