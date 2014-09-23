// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_filters/DdgFilter.cxxtest.hh
/// @brief  test for ddG filter
/// @author Rocco Moretti (rmoretti@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/simple_filters/DdgFilter.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/conformation/Residue.hh>

// Utility Headers
#include <utility/file/FileName.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.simple_filters.DdgFilters.cxxtest.hh");

// --------------- Test Class --------------- //

class DdgFilter : public CxxTest::TestSuite {

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
		sf->push_back( 212.0 ); // Bound no repack
		sf->push_back( 95 ); // Unbound no repack
		protocols::simple_filters::DdgFilter ddg_filter(0.0 /*threshold*/, scorefxn_, 1 /*jump*/, 1 /*repeats*/);

		(*scorefxn_)(*testpose_);
		ddg_filter.filter( sf );
		ddg_filter.repack_bound( false );
		ddg_filter.relax_bound( false );
    TS_ASSERT_EQUALS( ddg_filter.report_sm(*testpose_), 111 );

		ddg_filter.repack( false );
    TS_ASSERT_EQUALS( ddg_filter.report_sm(*testpose_), 117 );
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

		(*scorefxn_)(*testpose_);

    protocols::simple_filters::DdgFilter testfilter;
    TagCOP tag = tagptr_from_string("<Ddg name=test filter=sfT99 />\n");
    testfilter.parse_my_tag( tag, data, filters, movers, *testpose_ );

    TS_ASSERT_EQUALS( testfilter.report_sm(*testpose_), 99 );
  }
};
