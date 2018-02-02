// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/movers/LabelPoseFromResidueSelectorMover.cxxtest.hh
/// @brief test suite for protocols::movers::LabelPoseFromResidueSelectorMover
/// @author Jaume Bonet (jaume.bonet@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1rpb.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/fold_from_loops/movers/LabelPoseFromResidueSelectorMover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResiduePDBInfoHasLabelSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Utility Headers
#include <utility/file/FileName.hh>
#include <basic/Tracer.hh>

using namespace protocols::fold_from_loops::movers;

static basic::Tracer TR("protocols.fold_from_loops.movers.LabelPoseFromResidueSelectorMover.cxxtest.hh");

// --------------- Test Class --------------- //

class LabelPoseFromResidueSelectorMoverTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_name() {
		LabelPoseFromResidueSelectorMover lab;
		TS_ASSERT_EQUALS( lab.mover_name(), "LabelPoseFromResidueSelectorMover" );
	}

	void test_label() {
		LabelPoseFromResidueSelectorMover lab;
		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		core::select::residue_selector::ResidueIndexSelector selector( "1-5" );
		lab.residue_selector( selector );
		lab.label("TEST");
		lab.apply( trpcage );

		core::select::residue_selector::ResiduePDBInfoHasLabelSelector lsele("TEST");
		core::select::residue_selector::ResidueSubset subset = lsele.apply( trpcage );

		// test
		std::set < core::Size > acceptTrue;
		acceptTrue.insert( 1 );
		acceptTrue.insert( 2 );
		acceptTrue.insert( 3 );
		acceptTrue.insert( 4 );
		acceptTrue.insert( 5 );
		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT( !subset[ ii ] || ( acceptTrue.find( ii ) != acceptTrue.end() ) );
		}
	}

	void test_unlabel() {
		LabelPoseFromResidueSelectorMover lab1;
		LabelPoseFromResidueSelectorMover lab2;
		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		core::select::residue_selector::ResidueIndexSelector selector1( "1-5" );
		core::select::residue_selector::ResidueIndexSelector selector2( "3-5" );
		lab1.residue_selector( selector1 );
		lab1.label("TEST");
		lab1.apply( trpcage );
		lab2.residue_selector( selector2 );
		lab2.unlabel("TEST");
		lab2.apply( trpcage );

		core::select::residue_selector::ResiduePDBInfoHasLabelSelector lsele("TEST");
		core::select::residue_selector::ResidueSubset subset = lsele.apply( trpcage );

		// test
		std::set < core::Size > acceptTrue;
		acceptTrue.insert( 1 );
		acceptTrue.insert( 2 );
		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT( !subset[ ii ] || ( acceptTrue.find( ii ) != acceptTrue.end() ) );
		}
	}

};
