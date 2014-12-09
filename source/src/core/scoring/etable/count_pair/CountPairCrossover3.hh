// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/count_pair/CountPairCrossover3.hh
/// @brief  Count pair for residue pairs connected with one bond, where the
/// crossover from excluding to counting atom pair interactions is at 3 bonds.
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_scoring_etable_count_pair_CountPairCrossover3_hh
#define INCLUDED_core_scoring_etable_count_pair_CountPairCrossover3_hh

#include <core/scoring/etable/count_pair/CountPairFunction.hh>

#include <core/types.hh>


namespace core {
namespace scoring {
namespace etable {
namespace count_pair {

class CountPairCrossover3 : public CountPairFunction
{
public:
	public:
	typedef CountPairFunction parent;

public:

	virtual ~CountPairCrossover3();

	///@brief function used by derived classes and associated classes (like the trie's count pair hierarchy)
	static
	inline
	bool
	count_at_path_distance(
		Size const path_distance,
		Real & weight
	)
	{
		if ( path_distance < 3 ) return false;
		else if ( path_distance > 3 ) return true;

		weight = cp_half; return true;
	}

};

} // namespace count_pair
} // namespace etable
} // namespace scoring
} // namespace core

#endif
