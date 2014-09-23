// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_BaseBasePotential.fwd.hh
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_rna_RNA_LowResolutionPotential_fwd_hh
#define INCLUDED_core_scoring_rna_RNA_LowResolutionPotential_fwd_hh


// Unit headers
#include <utility/pointer/owning_ptr.hh>

// Package headers

// Project headers

// Utility headers

// C++


namespace core {
namespace scoring {
namespace rna {

class RNA_LowResolutionPotential;
typedef utility::pointer::shared_ptr< RNA_LowResolutionPotential > RNA_LowResolutionPotentialOP;

} //rna
} //scoring
} //core

#endif
