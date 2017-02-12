// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/free_peptide_movers/LongAxisRotationFreePeptideMover.cc
/// @brief Rotate the free peptide along the axis formed by its first and last CA atom
/// @author xingjiepan (xingjiepan@gmail.com)

#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/LongAxisRotationFreePeptideMover.hh>
#include <protocols/backbone_moves/local_backbone_mover/FreePeptide.hh>

#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.backbone_moves.local_backbone_mover.free_peptide_movers.LongAxisRotationFreePeptideMover" );


namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {
namespace free_peptide_movers {

LongAxisRotationFreePeptideMover::LongAxisRotationFreePeptideMover(Real radian, bool random):
 FreePeptideMover(), radian_(radian), random_(random) 
{
}

LongAxisRotationFreePeptideMover::~LongAxisRotationFreePeptideMover(){}

void
LongAxisRotationFreePeptideMover::apply(FreePeptide &free_peptide){
	using numeric::random::uniform;
	using numeric::rotation_matrix;

	Real rotation_radian = radian_;
	
	if(random_){
		rotation_radian = radian_ * uniform();	
		if(uniform() > 0.5){
			rotation_radian *= -1;
		}
	}

	xyzVector <Real> ca1 = free_peptide.ca_xyz(free_peptide.pivot1() + 1);
	xyzVector <Real> ca2 = free_peptide.ca_xyz(free_peptide.pivot2() - 1);

	free_peptide.rotate(rotation_matrix(ca2 - ca1, rotation_radian));
}


} //protocols
} //backbone_moves
} //local_backbone_mover
} //free_peptide_movers






