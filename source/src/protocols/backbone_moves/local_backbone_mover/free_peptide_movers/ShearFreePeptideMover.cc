// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/free_peptide_movers/ShearFreePeptideMover.cc
/// @brief Change the psi torsion of a residue and phi torsion of its subsequent residue by a opposite value.
/// @author xingjiepan (xingjiepan@gmail.com)

#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/ShearFreePeptideMover.hh>
#include <protocols/backbone_moves/local_backbone_mover/FreePeptide.hh>

#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/random/random.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.backbone_moves.local_backbone_mover.free_peptide_movers.ShearFreePeptideMover" );


namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {
namespace free_peptide_movers {

ShearFreePeptideMover::ShearFreePeptideMover(Size seqpos, Real torsion_degree, bool random):
	FreePeptideMover(), random_(random), seqpos_(seqpos), torsion_degree_(torsion_degree)
{

}

ShearFreePeptideMover::~ShearFreePeptideMover(){}

void 
ShearFreePeptideMover::apply(FreePeptide &free_peptide){
	using numeric::random::uniform;

	Real torsion_degree = torsion_degree_;

	if(random_){
		torsion_degree = torsion_degree_ * uniform();
		if(uniform() > 0.5){
			torsion_degree *= -1;
		}
	}

	free_peptide.psi(seqpos_, free_peptide.psi(seqpos_) + torsion_degree);
	free_peptide.phi(seqpos_ + 1, free_peptide.phi(seqpos_ + 1) - torsion_degree);
	free_peptide.align();
}

} //protocols
} //backbone_moves
} //local_backbone_mover
} //free_peptide_movers






