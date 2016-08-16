// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief general functions for generating typical kinds of Loops sets.
/// @author ashworth

#ifndef INCLUDED_protocols_loops_make_loops_hh
#define INCLUDED_protocols_loops_make_loops_hh

#include <protocols/loops/Loops.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace loops {

void loops_around_residues(
	Loops & loops,
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & residue_indices,
	core::Size gapsize = 6,
	core::Size extend = 2
);

} //namespace loops
} //namespace protocols

#endif
