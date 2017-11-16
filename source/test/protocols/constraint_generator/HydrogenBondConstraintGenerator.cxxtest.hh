// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/HydrogenBondConstraintGenerator.cxxtest.hh
/// @brief test suite for protocols::constraint_generator::HydrogenBondConstraintGenerator
/// @author Tom Linsky (tlinsky at uw dot edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1rpb.hh>

// Protocol headers
#include <protocols/constraint_generator/HydrogenBondConstraintGenerator.hh>

// Core headers
#include <core/conformation/Residue.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
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
using namespace core::scoring::func;
using namespace protocols::constraint_generator;

static basic::Tracer TR( "protocols.constraint_generator.HydrogenBondConstraintGenerator.cxxtest.hh" );

class HydrogenBondConstraintGeneratorTests : public CxxTest::TestSuite {

public:
	void setUp()
	{
		protocols_init();
	}


	void test_constraints_allresidues()
	{
		core::select::residue_selector::ResidueSelectorCOP selector1( new core::select::residue_selector::ResidueIndexSelector( "9" ) );
		core::select::residue_selector::ResidueSelectorCOP selector2( new core::select::residue_selector::ResidueIndexSelector( "16" ) );

		HydrogenBondConstraintGenerator hb_gen;
		hb_gen.set_id( "ap_generator1" );
		hb_gen.set_residue_selector1( selector1 );
		hb_gen.set_residue_selector2( selector2 );
		hb_gen.set_atom_pair_func( FuncOP( new FlatHarmonicFunc( 2.0, 0.5, 1.5 ) ) );
		hb_gen.set_angle1_func( FuncOP( new FlatHarmonicFunc( 2.09, 0.5, 0.4 ) ) );
		hb_gen.set_angle2_func( "FLAT_HARMONIC 1.90 0.5 0.4" );
		hb_gen.set_bounded( true );
		TS_ASSERT_EQUALS( hb_gen.class_name(), "HydrogenBondConstraintGenerator" );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::scoring::constraints::ConstraintCOPs const csts = hb_gen.apply( trpcage );

		for ( ConstraintCOPs::const_iterator c=csts.begin(); c!=csts.end(); ++c ) {
			// no constraint ptrs should be null
			TS_ASSERT( *c );
			// csts should all be AmbiguousConstraints
			TS_ASSERT( utility::pointer::dynamic_pointer_cast< AmbiguousConstraint const >( *c ) );
		}

		ConstraintSet cst_set;
		cst_set.add_constraints( csts );

		test::UTracer UT( "protocols/constraint_generator/HydrogenBondConstraintGenerator_default.cst" );
		cst_set.show_definition( UT, trpcage );
	}

	void test_constraints_hydrophobics()
	{
		// LEU 2 can't hydrogen bond to anything
		core::select::residue_selector::ResidueSelectorCOP selector1( new core::select::residue_selector::ResidueIndexSelector( "2" ) );
		core::select::residue_selector::ResidueSelectorCOP selector2( new core::select::residue_selector::ResidueIndexSelector( "16" ) );

		HydrogenBondConstraintGenerator hb_gen;
		hb_gen.set_id( "ap_generator1" );
		hb_gen.set_residue_selector1( selector1 );
		hb_gen.set_residue_selector2( selector2 );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::scoring::constraints::ConstraintCOPs const csts = hb_gen.apply( trpcage );
		TS_ASSERT_EQUALS( csts.size(), 0 );
	}

	void test_atoms()
	{
		core::select::residue_selector::ResidueSelectorCOP selector1( new core::select::residue_selector::ResidueIndexSelector( "9" ) );
		core::select::residue_selector::ResidueSelectorCOP selector2( new core::select::residue_selector::ResidueIndexSelector( "16" ) );

		HydrogenBondConstraintGenerator hb_gen;
		hb_gen.set_id( "ap_generator1" );
		hb_gen.set_residue_selector1( selector1 );
		hb_gen.set_residue_selector2( selector2 );
		hb_gen.set_atom_pair_func( FuncOP( new FlatHarmonicFunc( 2.0, 0.5, 1.5 ) ) );
		hb_gen.set_angle1_func( "FLAT_HARMONIC 2.09 0.5 0.4" );
		hb_gen.set_angle2_func( FuncOP( new FlatHarmonicFunc( 1.90, 0.5, 0.4 ) ) );
		hb_gen.set_atoms1( "OD2" );
		hb_gen.set_atoms2( "NE,NH2" );
		hb_gen.set_bounded( true );
		TS_ASSERT_EQUALS( hb_gen.class_name(), "HydrogenBondConstraintGenerator" );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::scoring::constraints::ConstraintCOPs const csts = hb_gen.apply( trpcage );

		for ( ConstraintCOPs::const_iterator c=csts.begin(); c!=csts.end(); ++c ) {
			// no constraint ptrs should be null
			TS_ASSERT( *c );
			// csts should all be AmbiguousConstraints
			TS_ASSERT( utility::pointer::dynamic_pointer_cast< AmbiguousConstraint const >( *c ) );
		}

		ConstraintSet cst_set;
		cst_set.add_constraints( csts );

		test::UTracer UT( "protocols/constraint_generator/HydrogenBondConstraintGenerator_9_OD1_16_NH2.cst" );
		cst_set.show_definition( UT, trpcage );
	}

	void test_complicated()
	{
		core::pose::Pose pose;
		core::io::pdb::build_pose_from_pdb_as_is( pose,
			"protocols/constraint_generator/HydrogenBondConstraintGenerator_complicated.pdb.gz" );

		core::select::residue_selector::ResidueSelectorCOP nuc( new core::select::residue_selector::ResidueIndexSelector( "18" ) );
		core::select::residue_selector::ResidueSelectorCOP his( new core::select::residue_selector::ResidueIndexSelector( "49" ) );
		core::select::residue_selector::ResidueSelectorCOP nuc_support( new core::select::residue_selector::ResidueIndexSelector( "53" ) );
		core::select::residue_selector::ResidueSelectorCOP acid( new core::select::residue_selector::ResidueIndexSelector( "27" ) );
		core::select::residue_selector::ResidueSelectorCOP acid_support( new core::select::residue_selector::ResidueIndexSelector( "45" ) );

		HydrogenBondConstraintGenerator hb_gen1;
		hb_gen1.set_id( "ap_generator_nuc_his" );
		hb_gen1.set_residue_selector1( nuc );
		hb_gen1.set_residue_selector2( his );
		hb_gen1.set_atoms2( "NE2" );
		hb_gen1.set_atom_pair_sd( 0.2 );
		hb_gen1.set_angle_sd( 0.2 );

		HydrogenBondConstraintGenerator hb_gen2;
		hb_gen2.set_id( "ap_generator_nuc_support" );
		hb_gen2.set_residue_selector1( nuc );
		hb_gen2.set_residue_selector2( nuc_support );

		HydrogenBondConstraintGenerator hb_gen3;
		hb_gen3.set_id( "ap_generator_his_acid" );
		hb_gen3.set_residue_selector1( his );
		hb_gen3.set_residue_selector2( acid );
		hb_gen3.set_atoms1( "ND1" );
		hb_gen3.set_atoms2( "OE1,OD1" );
		hb_gen3.set_atom_pair_sd( 0.2 );
		hb_gen3.set_angle_sd( 0.2 );

		HydrogenBondConstraintGenerator hb_gen4;
		hb_gen4.set_id( "ap_generator_acid_support" );
		hb_gen4.set_residue_selector1( acid );
		hb_gen4.set_residue_selector2( acid_support );
		hb_gen4.set_atoms1( "OE2,OD2" );
		hb_gen4.set_atom_pair_sd( 0.2 );
		hb_gen4.set_angle_sd( 0.2 );

		ConstraintSet cst_set;
		cst_set.add_constraints( hb_gen1.apply( pose ) );
		cst_set.add_constraints( hb_gen2.apply( pose ) );
		cst_set.add_constraints( hb_gen3.apply( pose ) );
		cst_set.add_constraints( hb_gen4.apply( pose ) );

		test::UTracer UT( "protocols/constraint_generator/HydrogenBondConstraintGenerator_complicated.cst" );
		cst_set.show_definition( UT, pose );
		cst_set.show_definition( TR.Debug, pose );

	}

};

