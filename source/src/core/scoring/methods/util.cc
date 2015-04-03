// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/methods/util.cc
///
/// @brief utility methods for scoring.
/// @author James Thompson

// Unit headers
#include <core/scoring/methods/util.hh>
#include <basic/options/option.hh>

#include <basic/Tracer.hh>
#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>

#include <utility/exit.hh>


// option key includes

#include <basic/options/keys/abinitio.OptionKeys.gen.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

static thread_local basic::Tracer tr( "core.scoring.methods" );

core::Real get_residue_weight_by_ss(
	const char ss
) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	core::Real rsd_wt = 1.0;
	if ( ss == 'H' ) {
		rsd_wt = option[ abinitio::rsd_wt_helix  ]();
	} else if ( ss == 'E' ) {
		rsd_wt = option[ abinitio::rsd_wt_strand ]();
	} else if ( ss == 'L' ) {
		rsd_wt = option[ abinitio::rsd_wt_loop   ]();
	} else {
		tr.Error << "Error: don't recognize secondary structure character '" <<  ss << "' " << std::endl;
		rsd_wt = option[ abinitio::rsd_wt_loop   ]();
	}

	return rsd_wt;
} // get_residue_wt_by_ss

bool residues_interact(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	core::Real const interaction_cutoff
) {

	Real const rsd1_reach( rsd1.nbr_radius() ),
		rsd2_reach( rsd2.nbr_radius() );
	Distance const intxn_dist( rsd1_reach + rsd2_reach + interaction_cutoff );
	DistanceSquared const intxn_dist2( intxn_dist * intxn_dist );
	DistanceSquared const nbr_dist2(
		rsd1.xyz( rsd1.nbr_atom() ).distance_squared( rsd2.xyz( rsd2.nbr_atom() ) )
	);

	return ( nbr_dist2 < intxn_dist2 );
}

bool atoms_interact(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	core::id::AtomID const & id1,
	core::id::AtomID const & id2,
	core::Real const interaction_cutoff
) {
debug_assert( id1.rsd() == rsd1.seqpos() );
debug_assert( id2.rsd() == rsd2.seqpos() );

	Distance const dist(
		rsd1.xyz( id1.atomno() ).distance( rsd2.xyz( id2.atomno() ) )
	);

	return ( dist < interaction_cutoff );
}

} // namespace methods
} // namespace scoring
} // namespace core
