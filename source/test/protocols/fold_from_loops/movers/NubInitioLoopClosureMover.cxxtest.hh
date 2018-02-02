// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/movers/NubInitioLoopClosureMover.cxxtest.hh
/// @brief test suite for protocols::movers::NubInitioLoopClosureMover
/// @author Jaume Bonet (jaume.bonet@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1rpb.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/variant_util.hh>

#include <protocols/fold_from_loops/movers/LabelPoseFromResidueSelectorMover.hh>
#include <protocols/fold_from_loops/movers/NubInitioLoopClosureMover.hh>
#include <protocols/fold_from_loops/filters/RmsdFromResidueSelectorFilter.hh>

#include <protocols/moves/Mover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragmentIO.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResiduePDBInfoHasLabelSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/util.hh>

// Utility Headers
#include <utility/file/FileName.hh>
#include <basic/Tracer.hh>

using namespace protocols::fold_from_loops::movers;
using namespace protocols::fold_from_loops::filters;


static basic::Tracer TR("protocols.fold_from_loops.movers.NubInitioLoopClosureMover.cxxtest.hh");

// --------------- Test Class --------------- //

class NubInitioLoopClosureMoverTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init_with_additional_options("-mute protocols.fold_from_loops");
	}

	void tearDown() {
	}

	void test_name() {
		NubInitioLoopClosureMover nic;
		TS_ASSERT_EQUALS( nic.mover_name(), "NubInitioLoopClosureMover" );
	}

	void test_not_labeled() {
		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		NubInitioLoopClosureMover nic;
		nic.apply( trpcage );
		TS_ASSERT_EQUALS( nic.get_last_move_status(), protocols::moves::FAIL_BAD_INPUT );
	}

	void test_labeled_closed() {

		RmsdFromResidueSelectorFilter rmsd;
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::pose::Pose refpose = create_trpcage_ideal_pose();
		core::fragment::FragSetOP frags = core::fragment::FragmentIO().read_data( "protocols/fold_from_loops/movers/frags.200.3mers" );
		rmsd.reference_pose( refpose );

		core::select::residue_selector::ResidueIndexSelector motif( "1-9,17-20" );
		core::select::residue_selector::ResidueIndexSelector tmplate( "10-16" );
		LabelPoseFromResidueSelectorMover lab;
		lab.residue_selector( motif );
		lab.label("MOTIF");
		lab.apply( trpcage );
		lab.label("HOTSPOT");
		lab.apply( trpcage );
		lab.residue_selector( tmplate );
		lab.label( "TEMPLATE" );
		lab.apply( trpcage );

		NubInitioLoopClosureMover nic;
		nic.fragments( frags );
		nic.apply( trpcage );
		TS_ASSERT_EQUALS( nic.get_last_move_status(), protocols::moves::MS_SUCCESS );

		core::select::residue_selector::ResiduePDBInfoHasLabelSelector slab("LOOPCLOSURE");
		TS_ASSERT_EQUALS( core::select::residue_selector::has_any_true_selection(slab.apply(trpcage)), 0 );
	}

	void test_labeled_open() {

		RmsdFromResidueSelectorFilter rmsd;
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::pose::Pose refpose = create_trpcage_ideal_pose();
		core::fragment::FragSetOP frags = core::fragment::FragmentIO().read_data( "protocols/fold_from_loops/movers/frags.200.3mers" );
		rmsd.reference_pose( refpose );

		core::select::residue_selector::ResidueIndexSelector motif( "1-9,17-20" );
		core::select::residue_selector::ResidueIndexSelector tmplate( "10-16" );
		LabelPoseFromResidueSelectorMover lab;
		lab.residue_selector( motif );
		lab.label("MOTIF");
		lab.apply( trpcage );
		lab.label("HOTSPOT");
		lab.apply( trpcage );
		lab.residue_selector( tmplate );
		lab.label( "TEMPLATE" );
		lab.apply( trpcage );

		core::pose::add_variant_type_to_pose_residue( trpcage, core::chemical::CUTPOINT_LOWER, 11 );
		core::pose::add_variant_type_to_pose_residue( trpcage, core::chemical::CUTPOINT_UPPER, 12 );

		NubInitioLoopClosureMover nic;
		nic.fragments( frags );
		nic.trust( true );
		nic.apply( trpcage );
		TS_ASSERT_EQUALS( nic.get_last_move_status(), protocols::moves::MS_SUCCESS );

		core::select::residue_selector::ResiduePDBInfoHasLabelSelector slab("LOOPCLOSURE");
		TS_ASSERT_LESS_THAN( 0, core::select::residue_selector::has_any_true_selection(slab.apply(trpcage)) );
	}

};
