// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/NubInitioMover.cxxtest.hh
/// @brief test suite for protocols::NubInitioMover
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
#include <protocols/fold_from_loops/filters/RmsdFromResidueSelectorFilter.hh>
#include <protocols/fold_from_loops/NubInitioMover.hh>
#include <protocols/fold_from_loops/utils/Nub.hh>

#include <protocols/constraint_generator/AtomPairConstraintGenerator.hh>
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

using namespace protocols::fold_from_loops;
using namespace protocols::fold_from_loops::movers;
using namespace protocols::fold_from_loops::filters;
using namespace protocols::fold_from_loops::utils;


static basic::Tracer TR("protocols.fold_from_loops.movers.NubInitioMover.cxxtest.hh");

// --------------- Test Class --------------- //

class NubInitioMoverTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init_with_additional_options("-mute protocols.fold_from_loops protocols.abinitio");
	}

	void tearDown() {
	}

	void test_name() {
		NubInitioMover ni;
		TS_ASSERT_EQUALS( ni.mover_name(), "NubInitioMover" );
	}

	void test_one_motif() {
		core::pose::Pose tmplateps = create_trpcage_ideal_pose();
		core::pose::Pose motifps   = create_trpcage_ideal_pose();
		core::select::residue_selector::ResidueIndexSelector insert_motif( "10-16" );

		protocols::constraint_generator::AtomPairConstraintGenerator pair_gen;
		pair_gen.set_id( "ap_generator1" );
		pair_gen.set_ca_only( true );
		pair_gen.set_max_distance( 25.0 );
		pair_gen.set_min_seq_sep( 6 );
		tmplateps.add_constraints( pair_gen.apply( tmplateps ) );

		NubInitioMover ni;
		NubOP nub( new Nub );

		nub->reference_pose( motifps.clone() );
		nub->selector( insert_motif.clone() );
		nub->small_fragments( core::fragment::FragmentIO().read_data( "protocols/fold_from_loops/movers/frags.200.3mers" ) );
		nub->large_fragments( core::fragment::FragmentIO().read_data( "protocols/fold_from_loops/movers/frags.200.9mers" ) );

		ni.template_selector( insert_motif.clone() );
		ni.rmsd_include_motif( true );
		ni.nub( nub );

		ni.apply( tmplateps );

		core::select::residue_selector::ResiduePDBInfoHasLabelSelector tmpsel("TEMPLATE");
		core::select::residue_selector::ResiduePDBInfoHasLabelSelector mtfsel("MOTIF");
		TS_ASSERT_EQUALS( core::select::residue_selector::count_selected(mtfsel.apply(tmplateps)), 7 );
		TS_ASSERT_EQUALS( core::select::residue_selector::count_selected(tmpsel.apply(tmplateps)), 13 );

	}

	void test_two_motifs() {
		core::pose::Pose tmplateps = create_trpcage_ideal_pose();
		core::pose::Pose motifps   = create_trpcage_ideal_pose();
		core::select::residue_selector::ResidueIndexSelector insert_motif( "1-5,18-20" );

		protocols::constraint_generator::AtomPairConstraintGenerator pair_gen;
		pair_gen.set_id( "ap_generator1" );
		pair_gen.set_ca_only( true );
		pair_gen.set_max_distance( 25.0 );
		pair_gen.set_min_seq_sep( 6 );
		tmplateps.add_constraints( pair_gen.apply( tmplateps ) );

		NubInitioMover ni;
		NubOP nub( new Nub );

		nub->reference_pose( motifps.clone() );
		nub->selector( insert_motif.clone() );
		nub->small_fragments( core::fragment::FragmentIO().read_data( "protocols/fold_from_loops/movers/frags.200.3mers" ) );
		nub->large_fragments( core::fragment::FragmentIO().read_data( "protocols/fold_from_loops/movers/frags.200.9mers" ) );

		ni.template_selector( insert_motif.clone() );
		ni.rmsd_include_motif( true );
		ni.nub( nub );

		ni.apply( tmplateps );

		core::select::residue_selector::ResiduePDBInfoHasLabelSelector tmpsel("TEMPLATE");
		core::select::residue_selector::ResiduePDBInfoHasLabelSelector mtfsel("MOTIF");
		TS_ASSERT_EQUALS( core::select::residue_selector::count_selected(mtfsel.apply(tmplateps)), 8 );
		TS_ASSERT_EQUALS( core::select::residue_selector::count_selected(tmpsel.apply(tmplateps)), 12 );

	}

};
