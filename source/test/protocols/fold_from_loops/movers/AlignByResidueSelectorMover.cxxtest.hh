// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/movers/AlignByResidueSelectorMover.cxxtest.hh
/// @brief test suite for protocols::movers::AlignByResidueSelectorMover
/// @author Jaume Bonet (jaume.bonet@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1rpb.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/fold_from_loops/filters/RmsdFromResidueSelectorFilter.hh>
#include <protocols/fold_from_loops/movers/AlignByResidueSelectorMover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Utility Headers
#include <utility/file/FileName.hh>
#include <basic/Tracer.hh>

using namespace protocols::fold_from_loops::movers;
using namespace protocols::fold_from_loops::filters;

static basic::Tracer TR("protocols.fold_from_loops.movers.AlignByResidueSelectorMover.cxxtest.hh");

// --------------- Test Class --------------- //

class AlignByResidueSelectorMoverTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_name() {
		AlignByResidueSelectorMover ali;
		TS_ASSERT_EQUALS( ali.mover_name(), "AlignByResidueSelectorMover" );
	}

	void test_all_true() {
		RmsdFromResidueSelectorFilter rmsd;
		AlignByResidueSelectorMover ali;
		core::pose::Pose trpcage1 = create_trpcage_ideal_pose();
		core::pose::Pose trpcage2 = create_trpcage_ideal_pose();

		rmsd.reference_pose( trpcage2 );
		rmsd.superimpose( false );
		ali.reference_pose( trpcage2 );
		ali.apply( trpcage1 );
		TS_ASSERT_EQUALS( rmsd.compute( trpcage1 ), 0.0 );
	}

	void test_shift_aligned() {
		core::select::residue_selector::ResidueIndexSelector selector1( "1-10" );
		core::select::residue_selector::ResidueIndexSelector selector2( "2-11" );

		RmsdFromResidueSelectorFilter rmsd;
		AlignByResidueSelectorMover ali;
		core::pose::Pose trpcage1 = create_trpcage_ideal_pose();
		core::pose::Pose trpcage2 = create_trpcage_ideal_pose();

		rmsd.reference_pose( trpcage2 );
		rmsd.superimpose( false );
		rmsd.reference_selector( selector1 );
		rmsd.query_selector( selector2 );
		ali.reference_pose( trpcage2 );
		ali.reference_selector( selector1 );
		ali.query_selector( selector2 );
		ali.apply( trpcage1 );

		TS_ASSERT_DELTA( rmsd.compute( trpcage1 ), 2.0217, 0.01 );
	}

};
