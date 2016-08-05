// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/constraint_generator/AtomPairConstraintGenerator.cxxtest.hh
/// @brief test suite for protocols::constraint_generator::AtomPairConstraintGenerator
/// @author Tom Linsky (tlinsky at uw dot edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1rpb.hh>

// Protocol headers
#include <protocols/constraint_generator/AddConstraints.hh>
#include <protocols/constraint_generator/AtomPairConstraintGenerator.hh>

// Core headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
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

static THREAD_LOCAL basic::Tracer TR( "protocols.constraint_generator.AtomPairConstraintGenerator.cxxtest.hh" );

class AtomPairConstraintGeneratorTests : public CxxTest::TestSuite {

public:

	void setUp()
	{
		protocols_init();
	}


	void test_constraints_allresidues()
	{
		AtomPairConstraintGenerator pair_gen;
		pair_gen.set_id( "ap_generator1" );
		pair_gen.set_ca_only( true );
		pair_gen.set_max_distance( 12.0 );
		pair_gen.set_min_seq_sep( 8 );
		TS_ASSERT_EQUALS( pair_gen.class_name(), "AtomPairConstraintGenerator" );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::scoring::constraints::ConstraintCOPs const csts = pair_gen.apply( trpcage );
		trpcage.dump_pdb( "test.pdb" );

		for ( ConstraintCOPs::const_iterator c=csts.begin(); c!=csts.end(); ++c ) {
			// no constraint ptrs should be null
			TS_ASSERT( *c );
			// csts should all be AtomPairConstraints
			TS_ASSERT( utility::pointer::dynamic_pointer_cast< AtomPairConstraint const >( *c ) );
		}

		ConstraintSet cst_set;
		cst_set.add_constraints( csts );

		test::UTracer UT( "protocols/constraint_generator/AtomPairConstraintGenerator_input_default.cst" );
		cst_set.show_definition( UT, trpcage );
	}

	void test_max_dist()
	{
		AtomPairConstraintGenerator pair_gen;
		pair_gen.set_id( "ap_generator1" );
		pair_gen.set_ca_only( true );
		pair_gen.set_max_distance( 8.0 );
		pair_gen.set_min_seq_sep( 3 );

		core::pose::Pose pose = create_trpcage_ideal_pose();
		ConstraintSet cst_set;
		cst_set.add_constraints( pair_gen.apply( pose ) );

		test::UTracer UT( "protocols/constraint_generator/AtomPairConstraintGenerator_input_distance.cst" );
		cst_set.show_definition( UT, pose );
	}

	void test_residue_selector()
	{
		core::select::residue_selector::ResidueIndexSelector selector( "1-10,12-20" );

		AtomPairConstraintGenerator pair_gen;
		pair_gen.set_id( "ap_generator1" );
		pair_gen.set_ca_only( true );
		pair_gen.set_max_distance( 12.0 );
		pair_gen.set_min_seq_sep( 8 );
		pair_gen.set_residue_selector( selector );

		core::pose::Pose pose = create_trpcage_ideal_pose();
		ConstraintSet cst_set;
		cst_set.add_constraints( pair_gen.apply( pose ) );

		test::UTracer UT( "protocols/constraint_generator/AtomPairConstraintGenerator_input_selector.cst" );
		cst_set.show_definition( UT, pose );
	}

	void test_all_atom()
	{
		AtomPairConstraintGenerator pair_gen;
		pair_gen.set_id( "ap_generator1" );
		pair_gen.set_ca_only( false );
		pair_gen.set_max_distance( 12.0 );
		pair_gen.set_min_seq_sep( 8 );

		core::pose::Pose pose = create_trpcage_ideal_pose();
		ConstraintSet cst_set;
		cst_set.add_constraints( pair_gen.apply( pose ) );

		test::UTracer UT( "protocols/constraint_generator/AtomPairConstraintGenerator_input_ca_cb.cst" );
		cst_set.show_definition( UT, pose );
	}

	void test_native()
	{
		core::pose::Pose pose( pdb1rpb_pose() );
		core::pose::PoseOP native( new core::pose::Pose );
		core::import_pose::pose_from_file( *native, "protocols/relax/AtomCoordinateCstMover_native.pdb", core::import_pose::PDB_file );

		AtomPairConstraintGenerator pair_gen;
		pair_gen.set_id( "ap_generator1" );
		pair_gen.set_ca_only( true );
		pair_gen.set_max_distance( 12.0 );
		pair_gen.set_min_seq_sep( 8 );
		pair_gen.set_reference_pose( native );

		ConstraintSet cst_set;
		cst_set.add_constraints( pair_gen.apply( pose ) );
		test::UTracer UT( "protocols/constraint_generator/AtomPairConstraintGenerator_native_default.cst" );
		cst_set.show_definition( UT, pose );
	}

};

