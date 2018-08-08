// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/gap_solution_pickers/RandomGapSolutionPicker.cc
/// @brief Pick a solution randomly.
/// @author xingjiepan (xingjiepan@gmail.com)

#include <protocols/backbone_moves/local_backbone_mover/gap_solution_pickers/RandomGapSolutionPicker.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>

// STD
#include <cstdlib>

static basic::Tracer TR( "protocols.backbone_moves.local_backbone_mover.gap_solution_pickers.RandomGapSolutionPicker" );


namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {
namespace gap_solution_pickers {

RandomGapSolutionPicker::RandomGapSolutionPicker():
	GapSolutionPicker()
{

}

RandomGapSolutionPicker::~RandomGapSolutionPicker()= default;

Size
RandomGapSolutionPicker::pick(core::pose::Pose const &, FreePeptide const &,
	vector1<vector1<Real> > const & pivot_torsions, Size const
) const {
	if ( 0 == pivot_torsions.size() ) { return 0; }
	return numeric::random::random_range( 1, pivot_torsions.size() );
}



} //protocols
} //backbone_moves
} //local_backbone_mover
} //gap_solution_pickers






