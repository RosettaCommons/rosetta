// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/free_peptide_movers/CircularPermuteFreePeptideMover.cc
/// @brief Circularly permute the torsions on the free peptide.
/// @author xingjiepan (xingjiepan@gmail.com)

#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/CircularPermuteFreePeptideMover.hh>
#include <protocols/backbone_moves/local_backbone_mover/FreePeptide.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.backbone_moves.local_backbone_mover.free_peptide_movers.CircularPermuteFreePeptideMover" );


namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {
namespace free_peptide_movers {

CircularPermuteFreePeptideMover::CircularPermuteFreePeptideMover(Size off_set, bool direction):
	FreePeptideMover(), off_set_(off_set), direction_(direction)
{}

CircularPermuteFreePeptideMover::~CircularPermuteFreePeptideMover()= default;

void
CircularPermuteFreePeptideMover::apply(FreePeptide &free_peptide){
	Size pivot1 = free_peptide.pivot1();
	Size pivot2 = free_peptide.pivot2();
	Size length = pivot2 - pivot1 - 1;

	vector1 < vector1 <Real> > torsions;

	// Store the current torsions

	for ( Size i = pivot1 + 1; i <= pivot2 - 1; ++i ) {
		vector1 <Real> torsion_one_res(3);
		torsion_one_res[1] = free_peptide.phi(i);
		torsion_one_res[2] = free_peptide.psi(i);
		torsion_one_res[3] = free_peptide.omega(i);

		torsions.push_back( torsion_one_res );
	}

	// Do circular permutation

	for ( Size i = pivot1 + 1; i <= pivot2 - 1; ++i ) {
		Size source_seqpos = (i + off_set_ - pivot1 - 1) % length + pivot1 + 1;
		if ( !direction_ ) {
			source_seqpos = (int(i) - off_set_ - pivot1 - 1) % length + pivot1 + 1;
		}

		free_peptide.phi(i, torsions[source_seqpos - pivot1][1]);
		free_peptide.psi(i, torsions[source_seqpos - pivot1][2]);
		free_peptide.omega(i, torsions[source_seqpos - pivot1][3]);
	}

	// Align the free peptide

	free_peptide.align();
}

} //protocols
} //backbone_moves
} //local_backbone_mover
} //free_peptide_movers






