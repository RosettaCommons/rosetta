// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/RNA_ScoringUtil.hh
/// @brief
/// @author Rhiju

#ifndef INCLUDED_core_scoring_rna_RNA_ScoringUtil_hh
#define INCLUDED_core_scoring_rna_RNA_ScoringUtil_hh

#include <core/types.hh>


namespace core {
namespace scoring {
namespace rna {

void
get_fade_correction(
	Real const z,
	Real const cutoff_lower,
	Real const cutoff_upper,
	Real const fade_zone,
	Real & fade_value,
	Real & fade_deriv );

} //rna
} //scoring
} //core

#endif
