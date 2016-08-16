// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/etable/count_pair/CountPairCrossover34.hh
/// @brief  Count pair for residue pairs where the crossover from excluding
/// to counting atom pair interactions is at 3 bonds, but where the weighting
/// is different for 3 bonds versus 4.
/// @author Jim Havranek (havranek@genetics.wustl.edu)


#ifndef INCLUDED_core_scoring_etable_count_pair_CountPairCrossover34_hh
#define INCLUDED_core_scoring_etable_count_pair_CountPairCrossover34_hh

#include <core/scoring/etable/count_pair/CountPairFunction.hh>

#include <core/types.hh>


namespace core {
namespace scoring {
namespace etable {
namespace count_pair {

class CountPairCrossover34 : public CountPairFunction
{
public:
public:
	typedef CountPairFunction parent;

public:

	virtual ~CountPairCrossover34();

	/// @brief function used by derived classes and associated classes (like the trie's count pair hierarchy)
	static
	inline
	bool
	count_at_path_distance(
		Size const path_distance,
		Real & weight
	)
	{

		// These weights are the values for
		// the amoeba permanent atomic multipole
		// potential.

		if ( path_distance > 4 ) {
			weight = 1.0;
			return true;
		} else if ( path_distance == 4 ) {
			weight = 0.8;
			return true;
		} else if ( path_distance == 3 ) {
			weight = 0.4;
			return true;
		} else {
			return false;
		}
	}
};

} // namespace count_pair
} // namespace etable
} // namespace scoring
} // namespace core

#endif
