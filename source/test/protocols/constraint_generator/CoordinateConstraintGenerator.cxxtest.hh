// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/ConstraintGeneratorFactory.cxxtest.hh
/// @brief test suite for protocols::constraint_generator::ConstraintGeneratorFactory
/// @author Tom Linsky (tlinsky at uw dot edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1rpb.hh>

// Protocol headers
#include <protocols/constraint_generator/AddConstraints.hh>
#include <protocols/constraint_generator/CoordinateConstraintGenerator.hh>

// Core headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Boost headers
#include <boost/assign.hpp>

// C++ headers

using namespace core::scoring::constraints;
using namespace protocols::constraint_generator;

static basic::Tracer TR( "protocols.constraint_generator.CoordinateConstraintGenerator.cxxtest.hh" );

class CoordinateConstraintGeneratorTests : public CxxTest::TestSuite {

public:

	void setUp()
	{
		protocols_init();
	}


	void test_constraints_allresidues()
	{
		CoordinateConstraintGenerator coord_cst;
		coord_cst.set_id( "coord_generator1" );
		coord_cst.set_sidechain( false );
		TS_ASSERT_EQUALS( coord_cst.class_name(), "CoordinateConstraintGenerator" );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::scoring::constraints::ConstraintCOPs const csts = coord_cst.apply( trpcage );

		// should be # csts equal to residues in the pose * 4 atoms per residue, plus one atom for OXT
		// total: 81 csts
		TS_ASSERT_EQUALS( csts.size(), 4*trpcage.size() + 1 );

		for ( core::scoring::constraints::ConstraintCOPs::const_iterator c=csts.begin(); c!=csts.end(); ++c ) {
			// no constraint ptrs should be null
			TS_ASSERT( *c );
			// csts should all be CoordinateConstraints
			TS_ASSERT( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::CoordinateConstraint const >( *c ) );
		}

		// With sidechains, should be #csts equal to number of atoms
		coord_cst.set_sidechain( true );
		core::scoring::constraints::ConstraintCOPs const csts2 = coord_cst.apply( trpcage );
		core::Size atom_count = 0;
		for ( core::Size resid=1; resid<=trpcage.size(); ++resid ) {
			atom_count += trpcage.residue( resid ).nheavyatoms();
		}
		TS_ASSERT_EQUALS( csts2.size(), atom_count );

		// with CA-only, should be #csts equal to # residues
		coord_cst.set_ca_only( true );
		core::scoring::constraints::ConstraintCOPs const csts3 = coord_cst.apply( trpcage );
		TS_ASSERT_EQUALS( csts3.size(), trpcage.size() );

		// if we add the constraints, energy should be zero
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::coordinate_constraint, 1.0 );

		core::Real const no_cst_score = (*scorefxn)( trpcage );
		TS_ASSERT_DELTA( no_cst_score, 0.0, 1e-4 );
		trpcage.add_constraints( csts2 );
		core::Real const cst_score = (*scorefxn)( trpcage );
		TS_ASSERT_DELTA( no_cst_score, cst_score, 1e-4 );

		// tweaking the pose should make score worse
		core::pose::Pose trpcage2 = trpcage;
		trpcage2.set_phi( 10, 80 );
		trpcage2.set_psi( 10, 80 );
		core::Real const moved_cst_score = (*scorefxn)( trpcage2 );
		TS_ASSERT_LESS_THAN( cst_score, moved_cst_score );
	}

	void test_constraints_subset()
	{
		core::select::residue_selector::ResidueSelectorCOP selector( new core::select::residue_selector::ResidueIndexSelector( "5-15,19" ) );
		CoordinateConstraintGenerator coord_cst;
		coord_cst.set_id( "coord_generator1" );
		coord_cst.set_sidechain( false );
		coord_cst.set_residue_selector( selector );
		TS_ASSERT_EQUALS( coord_cst.class_name(), "CoordinateConstraintGenerator" );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::scoring::constraints::ConstraintCOPs const csts = coord_cst.apply( trpcage );

		// should be # csts equal to residues in the selection * 4 atoms per residue
		// total: 48 csts
		TS_ASSERT_EQUALS( csts.size(), 48 );

		core::select::residue_selector::ResidueSubset const subset = selector->apply( trpcage );

		core::Size const root_resid = trpcage.fold_tree().root();

		utility::vector1< core::Size > cst_counts( trpcage.size(), 0 );
		for ( core::scoring::constraints::ConstraintCOPs::const_iterator c=csts.begin(); c!=csts.end(); ++c ) {
			// no constraint ptrs should be null
			TS_ASSERT( *c );
			// csts should all be CoordinateConstraints
			TS_ASSERT( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::CoordinateConstraint const >( *c ) );
			// all constrained residues should be in the subset
			utility::vector1< core::Size > const residues = (*c)->residues();
			for ( utility::vector1< core::Size >::const_iterator resid=residues.begin(); resid!=residues.end(); ++resid ) {
				TR.Debug << "Residue constrained: " << *resid << " selected: " << subset[ *resid ] << std::endl;
				if ( *resid != root_resid ) {
					TS_ASSERT( subset[ *resid ] );
				}
				cst_counts[ *resid ] += 1;
			}
		}

		// should be 4 constraints per residue (for each BB atom)
		// all csts should contain the root residue
		for ( core::Size resid=1; resid<=trpcage.size(); ++resid ) {
			if ( resid == root_resid ) {
				TS_ASSERT_EQUALS( cst_counts[ resid ], 48 );
			} else if ( subset[ resid ] ) {
				TS_ASSERT_EQUALS( cst_counts[ resid ], 4 );
			} else {
				TS_ASSERT_EQUALS( cst_counts[ resid ], 0 );
			}
		}
	}

	// --------------- Test Cases From AtomCoordinateCstMover --------------- //

	void test_input_harmonic_default() {
		core::pose::Pose pose( pdb1rpb_pose() );
		core::pose::addVirtualResAsRoot( pose );

		CoordinateConstraintGeneratorOP generator( new CoordinateConstraintGenerator );
		generator->set_sd( 0.5 );

		ConstraintSet cst_set;
		cst_set.add_constraints( generator->apply( pose ) );

		test::UTracer UT( "protocols/relax/AtomCoordinateCstMover_input_harmonic_default.cst" );
		cst_set.show_definition( UT, pose );
	}

	void test_input_bounded_default() {
		core::pose::Pose pose( pdb1rpb_pose() );
		core::pose::addVirtualResAsRoot( pose );

		CoordinateConstraintGeneratorOP generator( new CoordinateConstraintGenerator );
		generator->set_sd( 0.5 );
		generator->set_bounded( true );
		generator->set_bounded_width( 0.25 );

		ConstraintSet cst_set;
		cst_set.add_constraints( generator->apply( pose ) );

		test::UTracer UT( "protocols/relax/AtomCoordinateCstMover_input_bounded_default.cst" );
		cst_set.show_definition( UT, pose );
	}

	void test_native_harmonic_default() {
		core::pose::Pose pose( pdb1rpb_pose() );
		core::pose::addVirtualResAsRoot( pose );

		core::pose::PoseOP native( new core::pose::Pose );
		core::import_pose::pose_from_file( *native, "protocols/relax/AtomCoordinateCstMover_native.pdb", core::import_pose::PDB_file );
		pose.dump_pdb( "test.pdb" );

		CoordinateConstraintGeneratorOP generator( new CoordinateConstraintGenerator );
		generator->set_sd( 0.5 );
		generator->set_reference_pose( native );

		ConstraintSet cst_set;
		cst_set.add_constraints( generator->apply( pose ) );
		test::UTracer UT( "protocols/constraint_generator/CoordinateConstraintGenerator_native_harmonic_default.cst" );
		cst_set.show_definition( UT, pose );
	}

	void test_input_harmonic_sidechain() {
		core::pose::Pose pose( pdb1rpb_pose() );
		core::pose::addVirtualResAsRoot( pose );

		CoordinateConstraintGeneratorOP generator( new CoordinateConstraintGenerator );
		generator->set_sd( 0.75 );
		generator->set_sidechain( true );

		ConstraintSet cst_set;
		cst_set.add_constraints( generator->apply( pose ) );

		test::UTracer UT( "protocols/relax/AtomCoordinateCstMover_input_harmonic_sidechain.cst" );
		cst_set.show_definition( UT, pose );
	}

	void test_native_harmonic_sidechain() {
		core::pose::Pose pose( pdb1rpb_pose() );
		core::pose::addVirtualResAsRoot( pose );

		core::pose::PoseOP native( new core::pose::Pose );
		core::import_pose::pose_from_file( *native, "protocols/relax/AtomCoordinateCstMover_native.pdb", core::import_pose::PDB_file );

		CoordinateConstraintGeneratorOP generator( new CoordinateConstraintGenerator );
		generator->set_sd( 0.75 );
		generator->set_sidechain( true );
		generator->set_reference_pose( native );

		ConstraintSet cst_set;
		cst_set.add_constraints( generator->apply( pose ) );

		test::UTracer UT( "protocols/relax/AtomCoordinateCstMover_native_harmonic_sidechain.cst" );
		cst_set.show_definition( UT, pose );
	}

	void test_integration_native() {
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/relax/1a19.pdb", core::import_pose::PDB_file );

		core::pose::addVirtualResAsRoot( pose );

		core::pose::PoseOP native( new core::pose::Pose );
		core::import_pose::pose_from_file( *native, "protocols/relax/1a19_trunc.pdb", core::import_pose::PDB_file );

		CoordinateConstraintGeneratorOP generator( new CoordinateConstraintGenerator );
		generator->set_sd( 0.5 );
		generator->set_sidechain( true );
		generator->set_reference_pose( native );
		generator->set_ambiguous_hnq( true );

		ConstraintSet cst_set;
		cst_set.add_constraints( generator->apply( pose ) );

		test::UTracer UT("protocols/relax/AtomCoordinateCstMover_integration_native.cst");
		cst_set.show_definition( UT, pose );
	}

};
