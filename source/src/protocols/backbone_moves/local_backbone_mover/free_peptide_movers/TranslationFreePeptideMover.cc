// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/free_peptide_movers/TranslationFreePeptideMover.cc
/// @brief Translate the free peptide by a cartesian vector.
/// @author xingjiepan (xingjiepan@gmail.com)

#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/TranslationFreePeptideMover.hh>
#include <protocols/backbone_moves/local_backbone_mover/FreePeptide.hh>

#include <basic/Tracer.hh>

// Numeric headers

#include <numeric/random/random.hh>
#include <numeric/random/random_xyz.hh>


static basic::Tracer TR( "protocols.backbone_moves.local_backbone_mover.free_peptide_movers.TranslationFreePeptideMover" );


namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {
namespace free_peptide_movers {

TranslationFreePeptideMover::TranslationFreePeptideMover(xyzVector <Real> v_translate):
	FreePeptideMover(), random_(false), v_translate_(v_translate)
{
}

TranslationFreePeptideMover::TranslationFreePeptideMover(Real max_amplitude):
	FreePeptideMover(), random_(true), max_amplitude_(max_amplitude)
{
}

TranslationFreePeptideMover::~TranslationFreePeptideMover(){}

void
TranslationFreePeptideMover::apply(FreePeptide &free_peptide){
	using numeric::random::uniform;
	using numeric::random::random_normal;

	if ( random_ ) {
		v_translate_ = max_amplitude_ * uniform() * random_normal();
	}

	free_peptide.translate(v_translate_);
}

} //protocols
} //backbone_moves
} //local_backbone_mover
} //free_peptide_movers






