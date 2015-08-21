// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMBondLengthScore.cc
/// @brief  Molecular mechanics bond length score class
/// @author Frank DiMaio (based on Colin Smith's MMBondAngle potential)

// Unit headers
#include <core/scoring/mm/MMBondLengthScore.hh>
#include <core/scoring/mm/MMBondLengthLibrary.hh>

// Project headers
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/scoring/ScoringManager.hh>

// Utility header
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <iostream>
#include <string>
#include <map>

#include <utility/vector1.hh>
#include <utility/keys/Key2Tuple.hh>


namespace core {
namespace scoring {
namespace mm {

/// @details Auto-generated virtual destructor
MMBondLengthScore::~MMBondLengthScore() {}


MMBondLengthScore::MMBondLengthScore() :
	mm_bondlength_library_( scoring::ScoringManager::get_instance()->get_MMBondLengthLibrary() )
{ }

MMBondLengthScore::MMBondLengthScore( MMBondLengthLibrary const & mmtl ) :
	mm_bondlength_library_( mmtl )
{ }

/// @details dscore takes Kb, b0, and a distance;
/// returns an energy.
Real
MMBondLengthScore::score( Real Kb, Real b0, Real d ) const
{
	Real const del( d - b0 );
	return Kb * del * del;
}

/// @details Score take a set of 2 mm atom types in the form of an mm_bondlength_pair and a distance in Ang, and returns
/// the energy.
Real MMBondLengthScore::score( mm_bondlength_atom_pair mm_atomtype_set, Real d ) const {
	Real score = 0;

	// lookup
	mm_bondlength_library_citer_pair pair = mm_bondlength_library_.lookup(
		mm_atomtype_set.key1(), mm_atomtype_set.key2());

	// calc score
	for ( mm_bondlength_library_citer i = pair.first, e = pair.second; i != e; ++i ) {
		score += this->score( (i->second).key1(), (i->second).key2(), d );
	}

	return score;
}

/// @details dscore takes Kb, b0, and a distance;
/// returns an energy derivative.
Real
MMBondLengthScore::dscore( Real Kb, Real b0, Real d ) const {
	return 2 * Kb * (d-b0);
}

Real MMBondLengthScore::dscore( mm_bondlength_atom_pair mm_atomtype_set, Real d ) const {
	Real dscore_dang = 0;

	// lookup
	mm_bondlength_library_citer_pair pair = mm_bondlength_library_.lookup(
		mm_atomtype_set.key1(), mm_atomtype_set.key2() );

	// calc score
	for ( mm_bondlength_library_citer i = pair.first, e = pair.second; i != e; ++i ) {
		dscore_dang += dscore( (i->second).key1(), (i->second).key2(), d );
	}

	return dscore_dang;
}


} // namespace mm
} // namespace scoring
} // namespace core
