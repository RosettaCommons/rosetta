// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/P_AA_ABEGO3.cc
/// @brief  Amino acid probability given an ABEGO (ramachandran bins) triplet/sequence, arrays and functions
/// @author imv@uw.edu

// Unit headers
#include <core/scoring/P_AA_ABEGO3.hh>

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

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <fstream>
//MaximCode:
static THREAD_LOCAL basic::Tracer TR( "core.scoring.P_AA_ABEGO3" );

namespace core {
namespace scoring {

// 5 possible abego types (A, B, E, G, O).
#define ABEGO_COUNT 5

// The frequency/energy table contains all possible combinations of ABEGO at 3 positions
// and a central amino acid, which must be canonical (non-canonicals don't have sufficent frequency data)
#define KEY_COUNT (ABEGO_COUNT * ABEGO_COUNT * ABEGO_COUNT * core::chemical::AA::num_canonical_aas)

static_assert(core::chemical::AA::num_canonical_aas == 20, "The score file is predicated on existence of 20 canonical AA types");

enum ABEGO_index {
	ABEGO_index_INVALID=0,
	ABEGO_index_A,
	ABEGO_index_START = ABEGO_index_A,
	ABEGO_index_B,
	ABEGO_index_E,
	ABEGO_index_G,
	ABEGO_index_O,
	ABEGO_index_COUNT = ABEGO_index_O
};

inline ABEGO_index abego_to_index(char abego)
{
	switch(abego) {
	case 'A' : return ABEGO_index_A;
	case 'B' : return ABEGO_index_B;
	case 'E' : return ABEGO_index_E;
	case 'G' : return ABEGO_index_G;
	case 'O' : return ABEGO_index_O;
	};
	return ABEGO_index::ABEGO_index_INVALID;
}

/// @brief ctor -- Initialize the amino acid probability data structures
P_AA_ABEGO3::P_AA_ABEGO3() {
	read_P_AA_ABEGO3();
}

inline int abego3aa_to_index(const char abego1, const char abego2, const char abego3, const core::chemical::AA aa_index)
{
	ABEGO_index index1 = abego_to_index(abego1);
	ABEGO_index index2 = abego_to_index(abego2);
	ABEGO_index index3 = abego_to_index(abego3);

	if ( index1 == ABEGO_index_INVALID || index2 == ABEGO_index_INVALID || index3 == ABEGO_index_INVALID ) {
		TR << "Error: Invalid ABEGO triplet: '" << abego1 << abego2 << abego3 << std::endl;
		return -1;
	}

	if ( aa_index < core::chemical::first_l_aa || core::chemical::num_canonical_aas < aa_index ) {
		TR << "Error: Invalid AA, (int)aa='" << aa_index << "'" << std::endl;
		return -1;
	}

	// Base 5 encoding for abego section, * 20 canonical amino acids
	// Subtract 1 from each value because we're moving from 1-based enums
	int index = core::chemical::AA::num_canonical_aas * ( (index1-1) * 25 + (index2-1) * 5 + (index3-1) ) + (aa_index-1);
	return index;
}

inline int abego3aa_to_index(const char abego1, const char abego2, const char abego3, const char aa)
{
	TR.Debug << "looking up aa=" << aa << std::endl;
	core::chemical::AA aa_index = core::chemical::aa_from_oneletter_code(aa);
	TR.Debug << "aa_index returned " << aa << "=" << aa_index << std::endl;
	TR.Debug << "calling abego3aa_to_index(" << abego1 << "," << abego2 << "," << abego3 << "," << aa_index << ")" << std::endl;
	int index = abego3aa_to_index(abego1, abego2, abego3, aa_index);
	TR.Debug << "abego3aa_to_index returned " << index << std::endl;
	return index;
}

/// @brief Read the amino acid probability file into P_AA_ABEGO3
///
/// @note  Only the keys present in the file are given entries
void
P_AA_ABEGO3::read_P_AA_ABEGO3()
{
	TR << "loading ABEGO data file" << std::endl;
	//using namespace core::chemical;

	// Read the probability file and load the array

	p_AA_ABEGO3_.resize(KEY_COUNT);
	utility::io::izstream stream;
	basic::database::open( stream, "scoring/score_functions/p_aa_abego3.dat" );

	char aa_letter;
	char abego1;
	char abego2;
	char abego3;

	core::Real negative_energy; // Gabe's data is given as log(P1/P2) whereas energy is taken as negative, i.e. -log(P1/P2)

	while ( stream ) {
		// Each line lists 3 abego types, an amino acid letter, and the corresponding log probability, for example:
		// AGA N 0.499285392474
		// The ABEGO triplet is A, G, A, the amino acid in question is asparagine (corresponding to the G ABEGO),
		// its preference for this secondary structure are given as log( P(N|AGA)/P(N) ) = 0.499
		stream >> abego1 >> abego2 >> abego3 >> aa_letter >> negative_energy;

		int table_index = abego3aa_to_index(abego1, abego2, abego3, aa_letter);
		p_AA_ABEGO3_[table_index] = -negative_energy;
	}
	TR << "finished loading data file" << std::endl;

	stream.close();
}


/// @brief Returns -ln ( P(aa|abego3) / P(aa) ) where abego3 is the sequence of 3 ABEGO (ramachandran map) values,
/// centered on the aa in question. This measures how much the aa is conducive to (or can tolerate) being in
/// the current backbone conformation.
core::Real
P_AA_ABEGO3::P_AA_ABEGO3_energy( char abego_previous, char abego_current, char abego_next, core::chemical::AA aa ) const
{
	if ( chemical::num_canonical_aas < aa ) {
		return 0.0;
	}

	int index = abego3aa_to_index(abego_previous, abego_current, abego_next, aa);

	TR.Debug << "ABEGO3 index" << index << std::endl;
	core::Real ret = p_AA_ABEGO3_[index];
	return ret;
}


} // scoring
} // core

