// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/devel/denovo_design/CoreResiduesPerElementFilter.cxxtest.hh
/// @brief  test suite for devel::denovo_design::components::CoreResiduesPerElementFilter
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <devel/denovo_design/filters/CoreResiduesPerElementFilter.hh>

// Project headers

// Protocol headers
#include <protocols/matdes/SymDofMover.hh>
#include <protocols/moves/DsspMover.hh>

// Core headers
#include <core/kinematics/FoldTree.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// C++ headers
#include <protocols/denovo_design/test_utils.hh>

static basic::Tracer TR("devel.denovo_design.CoreResiduesPerElementFilter.cxxtest");

// --------------- Test Class --------------- //
class CoreResiduesPerElementFilterTests : public CxxTest::TestSuite {
	// scorefunction
	core::scoring::ScoreFunctionOP scorefxn;

public:

	// Shared initialization goes here.
	void setUp() {
		core_init();
		protocols_init();

		// set preserve header always
		basic::options::option[basic::options::OptionKeys::run::preserve_header].value(true);

		// initialize common filters/movers/scorefxns
		scorefxn = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction() );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_good_struct() {
		using namespace devel::denovo_design::filters;
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "devel/denovo_design/test_foldability.pdb" );

		CoreResiduesPerElementFilter corefilt;
		corefilt.set_core_cutoff( 2.0 );
		TS_ASSERT( corefilt.apply( input_pose ) );
	}

	void test_bad_struct() {
		using namespace devel::denovo_design::filters;
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/components/helix15.pdb" );

		CoreResiduesPerElementFilter corefilt;
		core::Real const badelements = corefilt.compute( input_pose );
		TS_ASSERT_EQUALS( badelements, 1 );
		TS_ASSERT( !corefilt.apply( input_pose ) );
	}

};
