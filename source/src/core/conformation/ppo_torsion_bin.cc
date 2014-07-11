// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/ppo_torsion_bin.hh
/// @brief  Functions for assigning and manipulating the phi/psi/omega torsion bins
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/conformation/ppo_torsion_bin.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

namespace core {
namespace conformation {

/// @details placed here to make it accessible to a wide range of applications, but it's quite possible
/// that placing this code elsewhere would be better bin boundaries are currently hard-coded -- ideally
/// in the future these can be read from external files, and thus adapted if desired
ppo_torsion_bin
get_torsion_bin (core::Real phi, core::Real  psi, core::Real omega) // omega defaults to 180 if not specified otherwise
{
	if (phi > 180) phi -= 360; // KIC returns positive values that need to be adjusted
	if (psi > 180) psi -= 360;
	if (omega > 180) omega -= 360;
	ppo_torsion_bin pos_bin = ppo_torbin_X;

	if (omega > 90 || omega < -90) { // trans --> uppercase letters
		if (phi <= 0) {
			if (psi < -130 || psi > 50) {
				pos_bin = ppo_torbin_B;
			} else { // -130 <= psi <= 50
				pos_bin = ppo_torbin_A;
			}
		} else { // phi > 0
			if (psi < -90 || psi > 90) {
				pos_bin = ppo_torbin_E;
			} else { // -90 <= psi <= 90
				pos_bin = ppo_torbin_G;
			}
		}
	} else { // cis --> lowercase letters
		if (phi <= 0) {
			if (psi < -130 || psi > 50) {
				pos_bin = ppo_torbin_b;
			} else { // -130 <= psi <= 50
				pos_bin = ppo_torbin_a;
			}
		} else { // phi > 0
			if (psi < -90 || psi > 90) {
				pos_bin = ppo_torbin_e;
			} else { // -90 <= psi <= 90
				pos_bin = ppo_torbin_g;
			}
		}
	}
	return pos_bin;
}

ppo_torsion_bin
remap_cis_omega_torsion_bins_to_trans( ppo_torsion_bin torbin ) {
	ppo_torsion_bin remapped_torbin;
	switch ( torbin ) {
		case ppo_torbin_A :
		case ppo_torbin_a :
			remapped_torbin = ppo_torbin_A;
			break;
		case ppo_torbin_B :
		case ppo_torbin_b :
			remapped_torbin = ppo_torbin_B;
			break;
		case ppo_torbin_E :
		case ppo_torbin_e :
			remapped_torbin = ppo_torbin_E;
			break;
		case ppo_torbin_G :
		case ppo_torbin_g :
			remapped_torbin = ppo_torbin_g;
			break;
		default :
			remapped_torbin = ppo_torbin_X;
	}
	return remapped_torbin;
}

bool
cis_omega_torsion_bin( ppo_torsion_bin torbin )
{
	switch (torbin) {
		case ppo_torbin_a:
		case ppo_torbin_b:
		case ppo_torbin_e:
		case ppo_torbin_g:
			return true;
		default:
			return false;
	}
}

torsion_bin_string
map_string_to_torsion_bin_string( std::string const & torstring )
{
	torsion_bin_string tbs( torstring.size() );
	for ( core::Size ii = 0; ii < torstring.size(); ++ii ) {
		if ( ! char_valid_as_torsion_bin( torstring[ ii ] )) {
			throw utility::excn::EXCN_Msg_Exception( "Invalid character '" + std::string(1,torstring[ii]) +
				"' at position " + utility::to_string( ii+1 ) + " of the torsion string \"" +
				torstring + "\" in trying to convert to a ppo_torsion_bin ");
		}
		tbs[ ii ] = map_char_to_torsion_bin( torstring[ ii ] );
	}
	return tbs;
}

bool
char_valid_as_torsion_bin( char torbin ) {
	switch ( torbin ) {
		case 'A' :
		case 'a' :
		case 'B' :
		case 'b' :
		case 'E' :
		case 'e' :
		case 'G' :
		case 'g' :
		case 'X' :
		case 'U' :
			return true;
		default  :
			return false;
	}
}

ppo_torsion_bin
map_char_to_torsion_bin( char torbin )
{
	switch ( torbin ) {
		case 'A' : return ppo_torbin_A;
		case 'a' : return ppo_torbin_a;
		case 'B' : return ppo_torbin_B;
		case 'b' : return ppo_torbin_b;
		case 'E' : return ppo_torbin_E;
		case 'e' : return ppo_torbin_e;
		case 'G' : return ppo_torbin_G;
		case 'g' : return ppo_torbin_g;
		case 'X' : return ppo_torbin_X;
		default  : return ppo_torbin_U;
	}
	// unreachable
	return ppo_torbin_U;
}

char
map_torsion_bin_to_char( ppo_torsion_bin torbin )
{
	switch ( torbin ) {
		case  ppo_torbin_A: return 'A';
		case  ppo_torbin_a: return 'a';
		case  ppo_torbin_B: return 'B';
		case  ppo_torbin_b: return 'b';
		case  ppo_torbin_E: return 'E';
		case  ppo_torbin_e: return 'e';
		case  ppo_torbin_G: return 'G';
		case  ppo_torbin_g: return 'g';
		case  ppo_torbin_X: return 'X';
		default : return 'U';
	}
	// unreachable
	return 'U';
}

} // namespace conformation
} // namespace core
