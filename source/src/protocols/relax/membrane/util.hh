// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/relax/membrane/util.hh
/// @brief      utility functions for the MPMutateRelaxMover
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_relax_membrane_util_hh
#define INCLUDED_protocols_relax_membrane_util_hh

// Unit Headers
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/scoring/ScoreFunction.fwd.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <basic/Tracer.fwd.hh>

namespace protocols {
namespace relax {
namespace membrane {

////////////////////////////////////////////////////////////////////////////////
/*
THIS IS HOW THE INPUT FORMAT OF THE FILE LOOKS LIKE:
= each line belongs to a single construct, i.e. a single sequence
= a single entry (format A163F) is a single mutation, multiple mutations per construct are possible
= example input:

A163F
W4N R27G K94E L45P
G32V P34N

= this means the first run is carried out for the single point mutation A163F
= the second run of the mover is carried out with a quadrupel mutation
= the third run is carried out with a double mutation
... and so on.

*/
////////////////////////////////////////////////////////////////////////////////

/// @brief Add mutants to private data: A163F into vectors
/// @details This is an entire line of the mutant input file
void add_mutant_to_vectors(
	core::pose::Pose & pose,
	std::string mutations,
	utility::vector1< utility::vector1< char > > & wt_res,
	utility::vector1< utility::vector1< core::Size > > & resid,
	utility::vector1< utility::vector1< char > > & mut_res );

/// @brief Check mutant file for errors
/// @details If Rosetta doesn't start crying, you're good to go
bool check_mutants_ok(
	core::pose::Pose & pose,
	utility::vector1< utility::vector1< char > > wt_res,
	utility::vector1< utility::vector1< core::Size > > resid );

} // membrane
} // relax
} // protocols

#endif // INCLUDED_protocols_relax_membrane_MPRangeRelaxMover_hh
