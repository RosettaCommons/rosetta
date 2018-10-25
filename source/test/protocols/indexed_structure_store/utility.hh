// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file test/core/indexed_structure_store/utility.hh
/// @brief Test utility functions for indexed_structure_store
/// @author Alex Ford (fordas@uw.edu)

#pragma once

#include <protocols/indexed_structure_store/Datatypes.hh>

#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/combine.hpp>

#include <cxxtest/TestSuite.h>

namespace protocols { namespace indexed_structure_store {


#define ASSERT_ALMOST_EQUAL(x, y, eps) TS_ASSERT_EQUALS(std::isnan(x), std::isnan(y)); if (!std::isnan(x) && !std::isnan(y)) { TS_ASSERT_DELTA(x, y, eps); }

inline void assert_residue_entry_almost_equal(
	ResidueEntry a, ResidueEntry b,
	bool check_structure_id = true,
	bool check_residue_id = true,
	bool check_sc = true,
	bool check_bb = true,
	bool check_orient = true,
	bool check_chain = true
) {
	if ( check_structure_id ) {
		TS_ASSERT_EQUALS(a.structure_id, b.structure_id);
	}

	if ( check_residue_id ) {
		TS_ASSERT_EQUALS(a.residue_id, b.residue_id);
	}


	if ( check_sc ) {
		ASSERT_ALMOST_EQUAL(a.sc.chi1, b.sc.chi1, 1e-3);
		ASSERT_ALMOST_EQUAL(a.sc.chi2, b.sc.chi2, 1e-3);
		ASSERT_ALMOST_EQUAL(a.sc.chi3, b.sc.chi3, 1e-3);
		ASSERT_ALMOST_EQUAL(a.sc.chi4, b.sc.chi4, 1e-3);
		TS_ASSERT_EQUALS(a.sc.aa, b.sc.aa);
	}

	if ( check_bb ) {
		ASSERT_ALMOST_EQUAL(a.bb.phi, a.bb.phi, 1e-3);
		ASSERT_ALMOST_EQUAL(a.bb.psi, a.bb.psi, 1e-3);
		ASSERT_ALMOST_EQUAL(a.bb.omega, a.bb.omega, 1e-3);
	}

	if ( check_orient ) {
		ASSERT_ALMOST_EQUAL(a.orient.N[0], b.orient.N[0], 1e-3);
		ASSERT_ALMOST_EQUAL(a.orient.N[1], b.orient.N[1], 1e-3);
		ASSERT_ALMOST_EQUAL(a.orient.N[2], b.orient.N[2], 1e-3);

		ASSERT_ALMOST_EQUAL(a.orient.C[0], b.orient.C[0], 1e-3);
		ASSERT_ALMOST_EQUAL(a.orient.C[1], b.orient.C[1], 1e-3);
		ASSERT_ALMOST_EQUAL(a.orient.C[2], b.orient.C[2], 1e-3);

		ASSERT_ALMOST_EQUAL(a.orient.CA[0], b.orient.CA[0], 1e-3);
		ASSERT_ALMOST_EQUAL(a.orient.CA[1], b.orient.CA[1], 1e-3);
		ASSERT_ALMOST_EQUAL(a.orient.CA[2], b.orient.CA[2], 1e-3);

		ASSERT_ALMOST_EQUAL(a.orient.O[0], b.orient.O[0], 1e-3);
		ASSERT_ALMOST_EQUAL(a.orient.O[1], b.orient.O[1], 1e-3);
		ASSERT_ALMOST_EQUAL(a.orient.O[2], b.orient.O[2], 1e-3);
	}

	if ( check_chain ) {
		TS_ASSERT_EQUALS(a.chain_ending, b.chain_ending);
	}
}

template< typename ResA, typename ResB>
void assert_residue_entries_almost_equal(
	ResA a, ResB b,
	bool check_structure_id = true,
	bool check_residue_id = true,
	bool check_sc = true,
	bool check_bb = true,
	bool check_orient = true,
	bool check_chain = true
) {
	TS_ASSERT_EQUALS( a.end() - a.begin(), b.end() - b.begin());

	for ( const auto & pair : boost::combine(a, b) ) {
		ResidueEntry ra, rb;
		boost::tie(ra, rb) = pair;
		assert_residue_entry_almost_equal(
			ra, rb,
			check_structure_id, check_residue_id, check_sc, check_bb, check_orient, check_chain
		);
	}
}

}}
