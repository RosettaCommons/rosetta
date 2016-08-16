// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/methods/util.cc
///
/// @brief utility methods for scoring.
/// @author James Thompson

// Unit headers
#include <core/scoring/methods/util.hh>
#include <basic/options/option.hh>

#include <basic/Tracer.hh>
#include <core/types.hh>

#include <core/scoring/Energies.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/scoring/PolymerBondedEnergyContainer.hh>

#include <utility/exit.hh>


// option key includes

#include <basic/options/keys/abinitio.OptionKeys.gen.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

static THREAD_LOCAL basic::Tracer tr( "core.scoring.methods" );

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

/// @brief Given two residues that may or may not be connected, determine which of the two, if any,
/// is the lower one and which is the upper.
/// @details Inputs are rsd1 and rsd2; outputs are rsd1_is_lo and rsd2_is_lo.  Both will be false if
/// the residues aren't conventionally connected (i.e. the C of one connected to the N of the other).
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
determine_lo_and_hi_residues(
	core::pose::Pose const &pose,
	core::Size const rsd1,
	core::Size const rsd2,
	bool &res1_is_lo,
	bool &res2_is_lo
) {
	res1_is_lo = (
		pose.residue(rsd1).has_upper_connect() &&
		pose.residue(rsd1).residue_connection_partner( pose.residue(rsd1).upper_connect().index() ) == rsd2 &&
		pose.residue(rsd2).has_lower_connect() &&
		pose.residue(rsd2).residue_connection_partner( pose.residue(rsd2).lower_connect().index() ) == rsd1
	);
	res2_is_lo = (
		pose.residue(rsd2).has_upper_connect() &&
		pose.residue(rsd2).residue_connection_partner( pose.residue(rsd2).upper_connect().index() ) == rsd1 &&
		pose.residue(rsd1).has_lower_connect() &&
		pose.residue(rsd1).residue_connection_partner( pose.residue(rsd1).lower_connect().index() ) == rsd2
	);
}

/// @brief Determines whether a long-range energies container exists in the pose energies object.  If not,
/// creates a new one and appends the score type to it, if necessary.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
create_long_range_energy_container(
	core::pose::Pose &pose,
	core::scoring::ScoreType const scoretype,
	core::scoring::methods::LongRangeEnergyType const lr_type
) {
	using namespace methods;

	// create LR energy container
	core::scoring::Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == 0 ) {
		create_new_lre_container = true;
	} else {
		LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		PolymerBondedEnergyContainerOP dec( utility::pointer::static_pointer_cast< core::scoring::PolymerBondedEnergyContainer > ( lrc ) );
		if ( !dec || !dec->is_valid( pose ) ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		utility::vector1< ScoreType > s_types;
		s_types.push_back( scoretype );
		LREnergyContainerOP new_dec( new PolymerBondedEnergyContainer( pose, s_types ) );
		energies.set_long_range_container( lr_type, new_dec );
	}
}

} // namespace methods
} // namespace scoring
} // namespace core
