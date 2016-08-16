// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/membrane/MPTerminiPenalty.fwd.hh
///
/// @brief  Membrane Protein Termini Penalty
/// @details Whole structure energy - penalty for residues on the wrong side of the membrane?
///    nd uses mpframework data
///    Last Modified: 3/31/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPTerminiPenalty_fwd_hh
#define INCLUDED_core_scoring_membrane_MPTerminiPenalty_fwd_hh

// Utility Methods
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace membrane {

class MPTerminiPenalty;
typedef utility::pointer::shared_ptr< MPTerminiPenalty > MPTerminiPenaltyOP;
typedef utility::pointer::shared_ptr< MPTerminiPenalty const > MPTerminiPenaltyCOP;

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_MPTerminiPenalty_fwd_hh
