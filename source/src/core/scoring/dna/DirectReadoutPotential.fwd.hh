// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/VDW_Energy.hh
/// @brief
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_dna_DirectReadoutPotential_fwd_hh
#define INCLUDED_core_scoring_dna_DirectReadoutPotential_fwd_hh

// Unit Headers

// Package headers

// Project headers

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace dna {

class DirectReadoutPotential;

typedef utility::pointer::shared_ptr< DirectReadoutPotential > DirectReadoutPotentialOP;
typedef utility::pointer::shared_ptr< DirectReadoutPotential const > DirectReadoutPotentialCOP;

}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
