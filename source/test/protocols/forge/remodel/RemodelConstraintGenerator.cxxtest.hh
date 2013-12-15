// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/forge/remodel/RemodelConstraintGenerator.cxxtest.hh
/// @brief  test suite for constraint generators
/// @author Tom Linsky

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <protocols/forge/constraints/RemoveCsts.hh>
#include <protocols/fldsgn/SheetConstraintsRCG.hh>
//#include <protocols/fldsgn/HSSTripletRCG.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/types.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>

//Auto Headers
#include <utility/vector1.hh>
using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("RemodelConstraintGenerator.cxxtest");

using namespace core;

class RemodelConstraintGeneratorTest : public CxxTest::TestSuite
{

public:
  RemodelConstraintGeneratorTest() {};

  core::scoring::ScoreFunctionOP scorefxn;
	protocols::fldsgn::SheetConstraintsRCGOP sheet_csts;
	protocols::forge::constraints::RemoveCsts rm_csts;
	protocols::fldsgn::SheetConstraintsRCGOP sheet_csts_badpair;
	//protocols::fldsgn::HSSTripletRCGOP hss_csts;

  // Shared initialization goes here.
  void setUp() {
    core_init();
    scorefxn = new scoring::ScoreFunction;
    scorefxn->reset();
	  scorefxn->set_weight( scoring::atom_pair_constraint, 1.0);
	  scorefxn->set_weight( scoring::angle_constraint, 1.0);
	  scorefxn->set_weight( scoring::dihedral_constraint, 1.0);
		sheet_csts = new protocols::fldsgn::SheetConstraintsRCG();
		sheet_csts->set_blueprint( "protocols/forge/remodel/test.blueprint" );
		sheet_csts_badpair = new protocols::fldsgn::SheetConstraintsRCG();
		sheet_csts_badpair->set_blueprint( "protocols/forge/remodel/bad_pair.blueprint" );
		TR << "Setting generator for rm_csts" << std::endl;
		rm_csts.set_generator( sheet_csts );
		TR << "Done setting generator for rm_csts" << std::endl;
		//readd_csts.cst_action( protocols::flxbb::ADD_PREGENERATED );
		//readd_csts.set_blueprint( "protocols/flxbb/test.blueprint" );
		//hss_csts = new protocols::fldsgn::HSSTripletRCG();
		//hss_csts->set_blueprint( "protocols/forge/remodel/test.blueprint" );
  }

  // Shared finalization goes here.
  void tearDown() {
  }

	/*// testing of HSSTripletRCG
	void Ntest_hsstriplet_constraints()
	{
		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, "protocols/forge/remodel/test.pdb" );
		core::Real score( scorefxn->score( pose ) );
		TS_ASSERT_DELTA( pose.energies().total_energies()[ scoring::atom_pair_constraint ], 0.0, 1e-4 );

		// now add constraints
		hss_csts->set_weight( 1.0 );
		hss_csts->set_max_distance( 12.0 );
		hss_csts->apply( pose );
		scorefxn->score( pose );
		scorefxn->show( TR, pose ); TR.flush();

		core::scoring::constraints::ConstraintCOPs csts( pose.constraint_set()->get_all_constraints() );
		// Should be two constraints for each helical residue -- helix sizes are 14,14,18,14
		TS_ASSERT_EQUALS( csts.size(), (core::Size)(14*2+14*2+18*2+14*2) );
		core::Size count( 0 );
		for ( core::Size i=1; i<=csts.size(); ++i ) {
			// type should be Ambiguous Constraint
			TS_ASSERT( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::AmbiguousConstraint const >( csts[i] ) );

			utility::vector1< core::Size > residues( csts[i]->residues() );
			for ( core::Size j=1; j<=residues.size(); ++j ) {
				if (residues[j] == 102) {
					++count;
				}
			}
		}
		// residue 102 should have 2 csts
		TS_ASSERT_EQUALS( count, (core::Size)2 );
		// total constraint energy should be 5.241
		TS_ASSERT_DELTA( pose.energies().total_energies()[ scoring::atom_pair_constraint ], 5.2408, 1e-4 );
	}*/

	// testing of SheetConstraintsRCG
  void test_sheet_constraints()
  {
    // first, we should try scoring a pose without constraints
    core::pose::Pose pose_nocst;
    core::import_pose::pose_from_pdb(pose_nocst, "protocols/forge/remodel/test.pdb" );
    //core::Real no_cst_score( scorefxn->score( pose_nocst ) );
		TR << "Before adding csts" << std::endl;
    //scorefxn->show( TR, pose_nocst );		TR.flush();
    core::Real dihedral_cst( pose_nocst.energies().total_energies()[ scoring::dihedral_constraint ] );
		TS_ASSERT_DELTA( dihedral_cst, 0.0, 1e-4 );

    core::pose::Pose testpose;
    core::import_pose::pose_from_pdb( testpose, "protocols/forge/remodel/test.pdb" );
    TS_ASSERT_EQUALS( testpose.constraint_set()->get_all_constraints().size(), (core::Size)0 );
		sheet_csts->set_angle_tolerance( 0.9 );
		sheet_csts->set_cacb_dihedral_tolerance( 0.9 );
		sheet_csts->set_bb_dihedral_tolerance( 0.9 );
    sheet_csts->apply( testpose );
		// check to make sure the proper number of constraints are there
		TR << "number of constraints after adding is " << testpose.constraint_set()->get_all_constraints().size() << std::endl;
		// There should be 31 constraints added to the pose
		// 16 atom pair constraints
		// 15 cacb dihedral constraints (not 16 because there is a GLY)
		// 32 bb dihedral constraints
		// 32 angle constraints
    TS_ASSERT_EQUALS( testpose.constraint_set()->get_all_constraints().size(), (core::Size)95 );
		// check for proper values of the constraints
    scorefxn->score( testpose );
    //scorefxn->show( TR, testpose );		TR.flush();

		// here, the constraint penalties should be 0 because the defaults are very forgiving
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::atom_pair_constraint ], 0.0, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::angle_constraint ], 0.0, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::dihedral_constraint ], 0.7491, 1e-4 );

		///////////////////////////////////////////////////////////////////////////////
		// now, import a new pose and use tighter constraints
		TS_ASSERT( testpose.remove_constraints() );
		testpose.clear();
		core::import_pose::pose_from_pdb( testpose, "protocols/forge/remodel/test.pdb" );
    TS_ASSERT_EQUALS( testpose.constraint_set()->get_all_constraints().size(), (core::Size)0 );

		// set constraints to tighter values than defaults
		sheet_csts->set_distance( 4.0 );
		sheet_csts->set_angle_tolerance( 0.35 );
		sheet_csts->set_cacb_dihedral_tolerance( 0.55 );
		sheet_csts->set_bb_dihedral_tolerance( 0.52 );
    sheet_csts->apply( testpose );
		// check to make sure the proper number of constraints are there
		TR << "number of constraints after adding is " << testpose.constraint_set()->get_all_constraints().size() << std::endl;
		// There should be 31 constraints added to the pose
		// 16 atom pair constraints
		// 15 dihedral constraints (not 16 because there is a GLY)
    TS_ASSERT_EQUALS( testpose.constraint_set()->get_all_constraints().size(), (core::Size)95 );
		// check for proper values of the constraints
    scorefxn->score( testpose );
    //scorefxn->show( TR, testpose );		TR.flush();

		// here, the constraint penalties should not be 0
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::atom_pair_constraint ], 10.8699, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::angle_constraint ], 0.0, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::dihedral_constraint ], 6.1114, 1e-4 );

		///////////////////////////////////////////////////////////////////////////////
		// we should be able to remove the constraints.
		rm_csts.set_generator( sheet_csts );
		rm_csts.apply( testpose );
		TS_ASSERT_EQUALS( testpose.constraint_set()->get_all_constraints().size(), (core::Size)0 );
		scorefxn->score( testpose );
		//scorefxn->show( TR, testpose );		TR.flush();

		// check to make sure all values for constraints are properly 0
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::atom_pair_constraint ], 0.0, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::angle_constraint ], 0.0, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::dihedral_constraint ], 0.0, 1e-4 );

		///////////////////////////////////////////////////////////////////////////////
		// we should be able to re-add the constraints
		sheet_csts->apply( testpose );
		TS_ASSERT_EQUALS( testpose.constraint_set()->get_all_constraints().size(), (core::Size)95 );
		TR << "AFTER re-adding constraints" << std::endl;
		scorefxn->score( testpose );
		//scorefxn->show( TR, testpose );		TR.flush();

		// check to make sure all values for constraints are properly 0
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::atom_pair_constraint ], 10.8699, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::angle_constraint ], 0.0, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::dihedral_constraint ], 6.1114, 1e-4 );

		///////////////////////////////////////////////////////////////////////////////
		// when we ramp up the coefficient, the constraints created should be higher.
		rm_csts.apply( testpose );
		TS_ASSERT_EQUALS( testpose.constraint_set()->get_all_constraints().size(), (core::Size)0 );
		scorefxn->score( testpose );
		//scorefxn->show( TR, testpose );		TR.flush();
		// check to make sure all values for constraints are properly 0
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::atom_pair_constraint ], 0.0, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::angle_constraint ], 0.0, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::dihedral_constraint ], 0.0, 1e-4 );

		sheet_csts->set_weight( 5.0 );
		sheet_csts->apply( testpose );
		TS_ASSERT_EQUALS( testpose.constraint_set()->get_all_constraints().size(), (core::Size)95 );
		scorefxn->score( testpose );
		//scorefxn->show( TR, testpose);		TR.flush();
		// check to make sure all values for constraints are properly 0
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::atom_pair_constraint ], 54.3497, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::angle_constraint ], 0.0, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::dihedral_constraint ], 30.5570, 1e-4 );
		// reset weight back to 1.0
		sheet_csts->set_weight( 1.0 );

		///////////////////////////////////////////////////////////////////////////////
		// clear csts and test only adding atom pair csts
		TS_ASSERT( testpose.remove_constraints() );
		testpose.clear();
		core::import_pose::pose_from_pdb( testpose, "protocols/forge/remodel/test.pdb" );
		sheet_csts->set_constrain_dist_only( true );
		sheet_csts->apply( testpose );
		TS_ASSERT_EQUALS( testpose.constraint_set()->get_all_constraints().size(), (core::Size)16 );
		scorefxn->score( testpose );
		//scorefxn->show( TR, testpose );		TR.flush();
		// check to make sure all values for constraints are properly 0
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::atom_pair_constraint ],10.8699, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::angle_constraint ], 0.0, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::dihedral_constraint ], 0.0, 1e-4 );
		// set to add angles/dihedrals again
		sheet_csts->set_constrain_dist_only( false );

		///////////////////////////////////////////////////////////////////////////////
		// TEST how well this function plays with existing constraints
		// reset the pose
		TS_ASSERT( testpose.remove_constraints() );
		testpose.clear();
		core::import_pose::pose_from_pdb( testpose, "protocols/forge/remodel/test.pdb" );

		// add a random constraint
		core::scoring::constraints::BoundFuncOP bound_func( new core::scoring::constraints::BoundFunc( 0.0, 4.0, 1.0, "random_constraint" ) );
		core::scoring::func::ScalarWeightedFuncOP func( new core::scoring::func::ScalarWeightedFunc( 1.0, bound_func ) );
		// CB on PHE3 to CB on ILE44 dist=4.9
		core::id::AtomID atom1( testpose.residue_type( 3 ).atom_index( "CB" ), 3 );
		core::id::AtomID atom2( testpose.residue_type( 44 ).atom_index( "CB" ), 44 );
		core::scoring::constraints::ConstraintOPs csts;
		core::scoring::constraints::AtomPairConstraintOP cst( new core::scoring::constraints::AtomPairConstraint( atom1, atom2, func ) );
		TS_ASSERT( cst );
		csts.push_back( cst );
		testpose.add_constraints( csts );
		TS_ASSERT_EQUALS( testpose.constraint_set()->get_all_constraints().size(), (core::Size)1 );
		scorefxn->score( testpose );
 		//scorefxn->show( TR, testpose );	TR.flush();
		// check to make sure all values for constraints are properly 0
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::atom_pair_constraint ], 0.6606, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::angle_constraint ], 0.0, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::dihedral_constraint ], 0.0, 1e-4 );

		// now add sheet csts
		sheet_csts->apply( testpose );
		TS_ASSERT_EQUALS( testpose.constraint_set()->get_all_constraints().size(), (core::Size)96 );
		scorefxn->score( testpose );
		//scorefxn->show( TR, testpose );		TR.flush();
		// check to make sure all values for constraints are properly 0
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::atom_pair_constraint ],11.5306, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::angle_constraint ], 0.0, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::dihedral_constraint ], 6.1114, 1e-4 );

		// now add another random cst
		// LYS18 CB <--> ARG93 CB  dist=11.8
		core::id::AtomID atom1_2( testpose.residue_type( 18 ).atom_index( "CB" ), 18 );
		core::id::AtomID atom2_2( testpose.residue_type( 93 ).atom_index( "CB" ), 93 );
		core::scoring::constraints::ConstraintOPs csts2;
		cst = new core::scoring::constraints::AtomPairConstraint( atom1_2, atom2_2, func );
		TS_ASSERT( cst );
		csts2.push_back( cst );
		testpose.add_constraints( csts2 );
		TS_ASSERT_EQUALS( testpose.constraint_set()->get_all_constraints().size(), (core::Size)97 );
		scorefxn->score( testpose );
 		//scorefxn->show( TR, testpose );		TR.flush();
		// check to make sure all values for constraints are properly 0
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::atom_pair_constraint ],19.1303, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::angle_constraint ], 0.0, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::dihedral_constraint ], 6.1114, 1e-4 );

		// now remove csts and there should only be 2
		rm_csts.apply( testpose );
		TS_ASSERT_EQUALS( testpose.constraint_set()->get_all_constraints().size(), (core::Size)2 );
		scorefxn->score( testpose );
 		//scorefxn->show( TR, testpose );		TR.flush();
		// check to make sure all values for constraints are properly 0
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::atom_pair_constraint ], 8.2603, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::angle_constraint ], 0.0, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::dihedral_constraint ], 0.0, 1e-4 );

		// now add constraints for sheets that contain a bad strand pairing (register shift off by one)
		sheet_csts->apply( testpose );
		TS_ASSERT_EQUALS( testpose.constraint_set()->get_all_constraints().size(), (core::Size)97 );

		// Constraints in centroid mode should still work and be removable, even if they were added in non-centroid mode
		TR << "Switching to centroid!" << std::endl;
		protocols::simple_moves::SwitchResidueTypeSetMover to_centroid_mover( core::chemical::CENTROID );
		to_centroid_mover.apply( testpose );

		sheet_csts_badpair->apply( testpose );
		TS_ASSERT_EQUALS( testpose.constraint_set()->get_all_constraints().size(), (core::Size)146 );
		rm_csts.apply( testpose );
		TS_ASSERT_EQUALS( testpose.constraint_set()->get_all_constraints().size(), (core::Size)91 );

		// check to make sure all values for constraints are properly set
		scorefxn->score( testpose );

		//scorefxn->show( TR, testpose );    TR.flush();
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::atom_pair_constraint ], 9.5225, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::angle_constraint ], 20.1780, 1e-4 );
		TS_ASSERT_DELTA( testpose.energies().total_energies()[ core::scoring::dihedral_constraint ], 230.7730, 1e-4 );

		// Remove constraints should be able to switch ids and remove the badpair_constraints
		rm_csts.set_generator( sheet_csts_badpair );
		rm_csts.apply( testpose );

		TS_ASSERT_EQUALS( testpose.constraint_set()->get_all_constraints().size(), (core::Size)2 );

		TR.flush();
  }

};
