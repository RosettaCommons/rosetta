// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/etable/count_pair/CountPairFactory.fwd.hh
/// @brief
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_scoring_etable_count_pair_CountPairFactory_hh
#define INCLUDED_core_scoring_etable_count_pair_CountPairFactory_hh

// Unit headers
#include <core/scoring/etable/count_pair/CountPairFactory.fwd.hh>

// Package headers
#include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>

// Project headers
#include <core/conformation/Residue.fwd.hh>

namespace core {
namespace scoring {
namespace etable {
namespace count_pair {

class CountPairFactory {

public:
	static
	CountPairFunctionOP
	create_count_pair_function(
		conformation::Residue const &,
		conformation::Residue const &,
		CPCrossoverBehavior
	);

	static
	void
	create_count_pair_function_and_invoke(
		conformation::Residue const &,
		conformation::Residue const &,
		CPCrossoverBehavior,
		Invoker & invoker
	);

	static
	CPResidueConnectionType
	determine_residue_connection(
		conformation::Residue const & res1,
		conformation::Residue const & res2
	);

	static
	CountPairFunctionOP
	create_intrares_count_pair_function(
		conformation::Residue const &,
		CPCrossoverBehavior
	);

};

class Invoker {
public:
	Invoker() {}
	virtual ~Invoker() {}

	virtual void invoke( CountPairFunction const & cp ) = 0;

};

} // namespace count_pair
} // namespace etable
} // namespace scoring
} // namespace core

#endif
