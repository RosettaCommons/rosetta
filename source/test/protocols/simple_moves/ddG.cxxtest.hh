// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/ddG.cxxtest.hh
/// @brief  test for ddG mover
/// @author Rocco Moretti (rmoretti@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/simple_moves/ddG.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/conformation/Residue.hh>

// Utility Headers
#include <utility/file/FileName.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.simple_moves.ddG.cxxtest.hh");

// --------------- Test Class --------------- //

class ddG : public CxxTest::TestSuite {

private:
	core::pose::PoseOP testpose_;
	core::scoring::ScoreFunctionOP scorefxn_;
public:

	void setUp() {
		core_init();

		testpose_ = create_2res_1ten_2res_trp_cage_poseop(); //dimer structure
		scorefxn_ = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
	}

	void tearDown() {
	}

	void test_use_filter() {
		StubMultiFilterOP sf( new StubMultiFilter( false ) );
		sf->push_back( 211.0 ); // Bound
		sf->push_back( 100 ); // Unbound
		protocols::simple_moves::ddG ddg_mover(scorefxn_);

		(*scorefxn_)(*testpose_);
		ddg_mover.filter( sf );
		ddg_mover.repeats( 1 );
		ddg_mover.repack_bound( false );
		ddg_mover.relax_bound( false );
		ddg_mover.apply( *testpose_ );
    TS_ASSERT_EQUALS( ddg_mover.sum_ddG(), 111 );
	}

  void test_filter_parsing() {
    basic::datacache::DataMap data;
    Filters_map filters;
    Movers_map movers;

		prime_Data( data );
		StubMultiFilterOP sf( new StubMultiFilter( false ) );
		sf->push_back( 199.0 ); // Bound
		sf->push_back( 100 ); // Unbound
    filters["sfT99"] = sf;

    protocols::simple_moves::ddG testmover;
    TagCOP tag = tagptr_from_string("<ddG name=test filter=sfT99 />\n");
    testmover.parse_my_tag( tag, data, filters, movers, *testpose_ );
		testmover.apply( *testpose_ );

    TS_ASSERT_EQUALS( testmover.sum_ddG(), 99 );
  }
};
