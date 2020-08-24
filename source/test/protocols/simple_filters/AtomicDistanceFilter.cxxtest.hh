// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/AtomicDistanceFilter.cxxtest.hh
/// @brief test suite for protocols::simple_filters::AtomicDistanceFilter
/// @author Jack Maguire, jackmaguire1444@gmail.com


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Package headers
#include <protocols/simple_filters/AtomicDistanceFilter.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/rosetta_scripts/XmlObjects.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <math.h>

#include <core/select/residue_selector/ResidueIndexSelector.hh>
//#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

using namespace protocols::simple_filters;

static basic::Tracer TR("protocols.simple_filters.AtomicDistanceFilter");


class AtomicDistanceFilterTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void test_AtomicDistanceFilter() {
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/simple_filters/interface_hbonds_input.pdb" );

		AtomicDistanceFilter filter;

		//ATOM      1  N   ALA A   1       5.619   0.314  -7.784  1.00  0.00           N
		filter.set_selector_for_res1(
			utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >( "1" )
		);
		filter.set_atomname1( "N" );

		//ATOM     39  CB  SER A   3       4.471   3.565  -3.815  1.00  0.00           C
		filter.set_selector_for_res2(
			utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >( "3" )
		);
		filter.set_atomname2( "CB" );

		constexpr core::Real dx = 5.619 - 4.471;
		constexpr core::Real dy = 0.314 - 3.565;
		constexpr core::Real dz = 7.784 - 3.815;

		constexpr core::Real dx_sq = dx*dx;
		constexpr core::Real dy_sq = dy*dy;
		constexpr core::Real dz_sq = dz*dz;

		core::Real const expected_distance = sqrt( dx_sq + dy_sq + dz_sq );

		TS_ASSERT_DELTA( filter.compute( pose ), expected_distance, 0.01 );

	}

};
