// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/constraint_generator/SegmentedAtomPairConstraintGenerator.cxxtest.hh
/// @brief test suite for protocols::constraint_generator::SegmentedAtomPairConstraintGenerator
/// @author Jaume Bonet (jaume.bonet@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1rpb.hh>

// Protocol headers
#include <protocols/constraint_generator/AddConstraints.hh>
#include <protocols/constraint_generator/AtomPairConstraintGenerator.hh>
#include <protocols/fold_from_loops/constraint_generator/SegmentedAtomPairConstraintGenerator.hh>

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
using namespace protocols::fold_from_loops::constraint_generator;

static basic::Tracer TR( "protocols.fold_from_loops.SegmentedAtomPairConstraintGenerator.cxxtest.hh" );

class SegmentedAtomPairConstraintGeneratorTests : public CxxTest::TestSuite {

public:

	void setUp()
	{
		protocols_init();
	}

	void test_constraints_defaults_allresidues()
	{
		SegmentedAtomPairConstraintGenerator spair_gen;
		TS_ASSERT_EQUALS( spair_gen.class_name(), "SegmentedAtomPairConstraintGenerator" );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::scoring::constraints::ConstraintCOPs const csts = spair_gen.apply( trpcage );

		for ( ConstraintCOPs::const_iterator c=csts.begin(); c!=csts.end(); ++c ) {
			// no constraint ptrs should be null
			TS_ASSERT( *c );
			// csts should all be AtomPairConstraints
			TS_ASSERT( utility::pointer::dynamic_pointer_cast< AtomPairConstraint const >( *c ) );
		}

		ConstraintSet cst_set;
		cst_set.add_constraints( csts );

		test::UTracer UT( "protocols/fold_from_loops/constraint_generator/SegmentedAtomPairConstraintGenerator_input_default.cst" );
		cst_set.show_definition( UT, trpcage );
	}

	void test_constraints_defaults_inner()
	{
		core::select::residue_selector::ResidueIndexSelector selector( "1-10,12-20" );

		SegmentedAtomPairConstraintGenerator spair_gen;
		spair_gen.set_residue_selector( selector );
		spair_gen.set_inner_min_seq_sep( 5 );
		spair_gen.set_do_outer( false );
		TS_ASSERT_EQUALS( spair_gen.class_name(), "SegmentedAtomPairConstraintGenerator" );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ConstraintSet cst_set;
		cst_set.add_constraints( spair_gen.apply( trpcage ) );

		test::UTracer UT( "protocols/fold_from_loops/constraint_generator/SegmentedAtomPairConstraintGenerator_input_inner.cst" );
		cst_set.show_definition( UT, trpcage );
	}

	void test_constraints_defaults_outer()
	{
		core::select::residue_selector::ResidueIndexSelector selector( "1-10,12-20" );

		SegmentedAtomPairConstraintGenerator spair_gen;
		spair_gen.set_residue_selector( selector );
		spair_gen.set_do_inner( false );
		TS_ASSERT_EQUALS( spair_gen.class_name(), "SegmentedAtomPairConstraintGenerator" );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ConstraintSet cst_set;
		cst_set.add_constraints( spair_gen.apply( trpcage ) );

		test::UTracer UT( "protocols/fold_from_loops/constraint_generator/SegmentedAtomPairConstraintGenerator_input_outer.cst" );
		cst_set.show_definition( UT, trpcage );
	}

};
