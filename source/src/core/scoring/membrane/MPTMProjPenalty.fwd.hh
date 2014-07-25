// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/membrane/MPTMProjPenalty.fwd.hh
///
///	@brief		Membrane Protein TM Proj Penalty
///	@details	Whole structure energy - Penalty for unreasonable tm-helix length compared to predicted
///				helix length (from topology) and uses mpframework data
///				Last Modified: 4/3/14
///
///	@author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPTMProjPenalty_fwd_hh
#define INCLUDED_core_scoring_membrane_MPTMProjPenalty_fwd_hh

// Utility Methods
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace membrane {
	
class MPTMProjPenalty;
typedef utility::pointer::owning_ptr< MPTMProjPenalty > MPTMProjPenaltyOP;
typedef utility::pointer::owning_ptr< MPTMProjPenalty const > MPTMProjPenaltyCOP;
	
} // membrane
} // scoring
} // core


#endif // INCLUDED_core_scoring_membrane_MPTMProjPenalty_fwd_hh
