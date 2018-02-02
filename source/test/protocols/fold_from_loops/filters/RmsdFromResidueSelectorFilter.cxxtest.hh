// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/filters/RmsdFromResidueSelectorFilter.cxxtest.hh
/// @brief test suite for protocols::filters::RmsdFromResidueSelectorFilter
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

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Utility Headers
#include <utility/file/FileName.hh>
#include <basic/Tracer.hh>

using namespace protocols::fold_from_loops::filters;

static basic::Tracer TR("protocols.fold_from_loops.filters.RmsdFromResidueSelector.cxxtest.hh");

// --------------- Test Class --------------- //

class RmsdFromResidueSelectorFilterTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_name() {
		RmsdFromResidueSelectorFilter rmsd;
		TS_ASSERT_EQUALS( rmsd.class_name(), "RmsdFromResidueSelectorFilter" );
	}

	void test_all_true() {
		RmsdFromResidueSelectorFilter rmsd;
		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		rmsd.reference_pose( trpcage );
		TS_ASSERT_EQUALS( rmsd.compute( trpcage ), 0.0 );
	}

	void test_selection_aligned() {
		core::select::residue_selector::ResidueIndexSelector selector1( "1-10" );
		core::select::residue_selector::ResidueIndexSelector selector2( "11-20" );

		RmsdFromResidueSelectorFilter rmsd;
		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		rmsd.reference_pose( trpcage );
		rmsd.reference_selector( selector1 );
		rmsd.query_selector( selector2 );
		TS_ASSERT_DELTA( rmsd.compute( trpcage ), 4.2348, 0.01 );
	}

	void test_selection_non_aligned() {
		core::select::residue_selector::ResidueIndexSelector selector1( "1-10" );
		core::select::residue_selector::ResidueIndexSelector selector2( "11-20" );

		RmsdFromResidueSelectorFilter rmsd;
		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		rmsd.reference_pose( trpcage );
		rmsd.reference_selector( selector1 );
		rmsd.query_selector( selector2 );
		rmsd.superimpose( false );
		TS_ASSERT_DELTA( rmsd.compute( trpcage ), 13.0370, 0.01 );
	}

};
