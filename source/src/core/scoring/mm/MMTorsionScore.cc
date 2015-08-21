// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMTorsionScore.cc
/// @brief  Molecular mechanics torsion score class
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit headers
#include <core/scoring/mm/MMTorsionScore.fwd.hh>
#include <core/scoring/mm/MMTorsionScore.hh>
#include <core/scoring/mm/MMTorsionLibrary.fwd.hh>
#include <core/scoring/mm/MMTorsionLibrary.hh>

// Project headers
#include <core/chemical/MMAtomTypeSet.fwd.hh>

#include <core/scoring/ScoringManager.hh>

// Utility header
#include <utility/keys/Key4Tuple.hh>
#include <utility/keys/Key3Tuple.hh>

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <iostream>
#include <string>
#include <map>
#include <math.h>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace mm {

/// @details Auto-generated virtual destructor
MMTorsionScore::~MMTorsionScore() {}

MMTorsionScore::MMTorsionScore() :
	mm_torsion_library_( scoring::ScoringManager::get_instance()->get_MMTorsionLibrary() )
{ }

MMTorsionScore::MMTorsionScore( MMTorsionLibrary const & mmtl ) :
	mm_torsion_library_( mmtl )
{ }

/// @details Score take a set of 4 mm atom types in the form of an mm_torsion_quad and an angle in radians and returns
/// an energy. This is done by first looking up the set(s) of torsional parameters in the library and then calculating
/// the score given an angle in radians.
Real
MMTorsionScore::score( mm_torsion_atom_quad mm_atomtype_set, Real angle ) const
{
	Real score = 0;

	// lookup
	mm_torsion_library_citer_pair pair = mm_torsion_library_.lookup(
		mm_atomtype_set.key1(),
		mm_atomtype_set.key2(),
		mm_atomtype_set.key3(),
		mm_atomtype_set.key4() );

	// calc score
	for ( mm_torsion_library_citer i = pair.first, e = pair.second; i != e; ++i ) {
		score += ( (i->second).key1() * ( 1+cos( (i->second).key2() *  angle - (i->second).key3() ) ) );
	}

	/* Debug virtual atom scores
	if ( mm_atomtype_set.key1() == 38
	|| mm_atomtype_set.key2() == 38
	|| mm_atomtype_set.key3() == 38
	|| mm_atomtype_set.key4() == 38 ) {
	std::cout << "MM virtual score: " << score << std::endl;
	if ( score != 0.0 ) {
	for ( mm_torsion_library_citer i = pair.first, e = pair.second; i != e; ++i ) {
	std::cout << "(i->second).key1() " << (i->second).key1();
	std::cout << " (i->second).key2() " << (i->second).key2();
	std::cout << " (i->second).key3() " << (i->second).key3() << std::endl;
	}

	}
	}
	*/

	return score;
}

/// @details dscore take a set of 4 mm atom types in the form of an mm_torsion_quad and an angle in radians and returns
/// an energy derivative. This is done by first looking up the set(s) of torsional parameters in the library and then calculating
/// the score derivative given an angle in radians.
Real
MMTorsionScore::dscore( mm_torsion_atom_quad mm_atomtype_set, Real angle ) const
{
	Real dscore_dang = 0;

	// lookup
	mm_torsion_library_citer_pair pair = mm_torsion_library_.lookup(
		mm_atomtype_set.key1(),
		mm_atomtype_set.key2(),
		mm_atomtype_set.key3(),
		mm_atomtype_set.key4() );

	// calc score
	//for ( mm_torsion_library_citer i = pair.first, e = pair.second; i != e; ++i ) {
	// score += ( (i->second).key1() * ( 1+cos( (i->second).key2() *  angle - (i->second).key3() ) ) );
	//}

	for ( mm_torsion_library_citer i = pair.first, e = pair.second; i != e; ++i ) {
		/// The math below is entirely Doug's and I'm trusting it.
		dscore_dang += (-1 * (i->second).key1() * (i->second).key2() * sin( (i->second).key2() * angle - (i->second).key3() ) );
	}

	/* Debug virtual atom scores
	if ( mm_atomtype_set.key1() == 38
	|| mm_atomtype_set.key2() == 38
	|| mm_atomtype_set.key3() == 38
	|| mm_atomtype_set.key4() == 38 ) {
	std::cout << "MM virtual score: " << dscore_dang << std::endl;
	if ( score != 0.0 ) {
	for ( mm_torsion_library_citer i = pair.first, e = pair.second; i != e; ++i ) {
	std::cout << "(i->second).key1() " << (i->second).key1();
	std::cout << " (i->second).key2() " << (i->second).key2();
	std::cout << " (i->second).key3() " << (i->second).key3() << std::endl;
	}

	}
	}
	*/

	return dscore_dang;
}

} // namespace mm
} // namespace scoring
} // namespace core
