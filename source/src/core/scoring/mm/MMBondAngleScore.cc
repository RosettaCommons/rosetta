// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMBondAngleScore.cc
/// @brief  Molecular mechanics bond angle score class
/// @author Colin A. Smith (colin.smith@ucsf.edu)

// Unit headers
#include <core/scoring/mm/MMBondAngleScore.fwd.hh>
#include <core/scoring/mm/MMBondAngleScore.hh>
#include <core/scoring/mm/MMBondAngleLibrary.fwd.hh>
#include <core/scoring/mm/MMBondAngleLibrary.hh>

// Project headers
#include <core/chemical/MMAtomTypeSet.fwd.hh>

#include <core/scoring/ScoringManager.hh>

// Utility header
#include <utility/keys/Key3Tuple.hh>
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
MMBondAngleScore::~MMBondAngleScore() {}

MMBondAngleScore::MMBondAngleScore() :
	mm_bondangle_library_( scoring::ScoringManager::get_instance()->get_MMBondAngleLibrary() )
{ }

MMBondAngleScore::MMBondAngleScore( MMBondAngleLibrary const & mmtl ) :
	mm_bondangle_library_( mmtl )
{ }

/// @details score takes Ktheta, theta0, and an angle in radians and returns
/// an energy.
Real
MMBondAngleScore::score( Real Ktheta, Real theta0, Real angle ) const
{
	Real const angle_diff( angle - theta0 );
	return Ktheta * angle_diff * angle_diff;
}

/// @details Score take a set of 3 mm atom types in the form of an mm_bondangle_tri and an angle in radians and returns
/// an energy. This is done by first looking up the set(s) of bond angle parameters in the library and then calculating
/// the score given an angle in radians.
Real
MMBondAngleScore::score( mm_bondangle_atom_tri mm_atomtype_set, Real angle ) const
{
	Real score = 0;

	// lookup
	mm_bondangle_library_citer_pair pair = mm_bondangle_library_.lookup(
		mm_atomtype_set.key1(),
		mm_atomtype_set.key2(),
		mm_atomtype_set.key3() );

	// calc score
	for ( mm_bondangle_library_citer i = pair.first, e = pair.second; i != e; ++i ) {

		score += this->score( (i->second).key1(), (i->second).key2(), angle );

		//std::cout << "[" << mm_atomtype_set.key1() << "," <<mm_atomtype_set.key2() << "," << mm_atomtype_set.key3() << ": " << numeric::conversions::degrees( angle ) << " " << numeric::conversions::degrees((i->second).key2()) << " " <<  (i->second).key1() << ": " << (i->second).key1() * angle_diff * angle_diff << "]" << std::endl;
	}

	/* Debug virtual atom scores
	if ( mm_atomtype_set.key1() == 38
	|| mm_atomtype_set.key2() == 38
	|| mm_atomtype_set.key3() == 38 ) {
	std::cout << "MM virtual score: " << score << std::endl;
	if ( score != 0.0 ) {
	for ( mm_bondangle_library_citer i = pair.first, e = pair.second; i != e; ++i ) {
	std::cout << "(i->second).key1() " << (i->second).key1();
	std::cout << " (i->second).key2() " << (i->second).key2() << std::endl;
	}

	}
	}
	*/
	return score;
}

/// @details dscore takes Ktheta, theta0, and an angle in radians and returns
/// an energy derivative.
Real
MMBondAngleScore::dscore( Real Ktheta, Real theta0, Real angle ) const
{
	Real const angle_diff( angle - theta0 );
	return 2 * Ktheta * angle_diff;
}

/// @details dScore take a set of 3 mm atom types in the form of an mm_bondangle_tri and an angle in radians and returns
/// an energy derivative. This is done by first looking up the set(s) of bond angle parameters in the library and then calculating
/// the dscore_dang.
Real
MMBondAngleScore::dscore( mm_bondangle_atom_tri mm_atomtype_set, Real angle ) const
{
	Real dscore_dang = 0;

	// lookup
	mm_bondangle_library_citer_pair pair = mm_bondangle_library_.lookup(
		mm_atomtype_set.key1(),
		mm_atomtype_set.key2(),
		mm_atomtype_set.key3() );

	// calc score
	for ( mm_bondangle_library_citer i = pair.first, e = pair.second; i != e; ++i ) {
		dscore_dang += dscore( (i->second).key1(), (i->second).key2(), angle );
		//std::cout << "core.mm.MMBondAngleEnergy:     dscore_dang = " << dscore_dang <<
		//" sc: " << score( mm_atomtype_set, angle ) << " +d " << score( mm_atomtype_set, angle + 1e-4 ) <<
		//" -d " << score( mm_atomtype_set, angle - 1e-4 ) << " numderiv: " << (score( mm_atomtype_set, angle + 1e-4 ) - score( mm_atomtype_set, angle - 1e-4 ) ) / 2e-4 << std::endl;
	}

	return dscore_dang;
}

} // namespace mm
} // namespace scoring
} // namespace core
