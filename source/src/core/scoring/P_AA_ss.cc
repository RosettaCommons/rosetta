// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/P_AA_ss.cc
/// @brief  Amino acid probability arrays and functions
/// @author FD

// Unit headers
#include <core/scoring/P_AA_ss.hh>

// Project headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <core/conformation/Residue.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// Utility headers
#include <utility/io/izstream.hh>

// C++ headers
#include <utility/assert.hh>

#include <utility/vector1.hh>

//MaximCode:
static basic::Tracer TR( "core.scoring.P_AA_ss" );


namespace core {
namespace scoring {


/// @brief ctor -- Initialize the amino acid probability data structures
P_AA_ss::P_AA_ss() {
	read_P_AA_ss();
}


/// @brief Read the amino acid probability file into P_AA_ss
///
/// @note  Only the keys present in the file are given entries
void
P_AA_ss::read_P_AA_ss()
{
	using namespace core::chemical;

	// Read the probability file and load the array
	utility::io::izstream stream;
	basic::database::open( stream, "scoring/score_functions/p_aa_ss.dat" );  //fpd to do: make this a flag

	p_L_.resize( num_canonical_aas );
	p_H_.resize( num_canonical_aas );
	p_E_.resize( num_canonical_aas );
	p0_L_ = p0_H_ = p0_E_ = 0.0;

	while ( stream ) {
		std::string id;
		core::Real L_i, H_i, E_i;
		stream >> id >> L_i >> H_i >> E_i;
		if ( stream ) {
			if ( id == "XXX" ) {
				p0_L_ = L_i;
				p0_H_ = H_i;
				p0_E_ = E_i;
			} else {
				AA aa = aa_from_name( id );
				p_L_[ aa ] = L_i;
				p_H_[ aa ] = H_i;
				p_E_[ aa ] = E_i;
			}
		}
	}
	stream.close();
}


/// @brief Probability energies from P(aa|phi,psi): Low level calculation for non-terminus position
/// You must pass an L amino acid 1-20 to this function! - if res have a backbone aa, apply it first
/// and if res is D, switch to L first!
core::Real
P_AA_ss::P_AA_ss_energy( chemical::AA aa, char ss ) const
{
	if ( aa > chemical::num_canonical_aas ) return 0.0;

	if ( ss == 'L' ) return p_L_[aa] + p0_L_;
	if ( ss == 'H' ) return p_H_[aa] + p0_H_;
	if ( ss == 'E' ) return p_E_[aa] + p0_E_;

	TR << "Warning!  Unrecognized ss type " << ss << std::endl;
	return 0.0;
}


} // namespace scoring
} // namespace rosetta

