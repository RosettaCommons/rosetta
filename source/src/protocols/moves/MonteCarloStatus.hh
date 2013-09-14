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
/// @author

#ifndef INCLUDED_protocols_moves_MonteCarloStatus_hh
#define INCLUDED_protocols_moves_MonteCarloStatus_hh

#include <string>
namespace protocols {
namespace moves {

// 	mc_accepted
// 		3 = accepted:score beat low score and last_accepted score
// 		2 = accepted:score beat last_accepted score
// 		1 = thermally accepted: score worse than last_accepted score
// 		0 = not accepted
typedef enum {
	MCA_accepted_score_beat_low=3,
	MCA_accepted_score_beat_last=2,
	MCA_accepted_thermally=1,
	MCA_rejected=0
} MCA; // mc_accepted state

std::string
to_string( MCA const & mc_accepted );

} // simple_moves
} // protocols


#endif
