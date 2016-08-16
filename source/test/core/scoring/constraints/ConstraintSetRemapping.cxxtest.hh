// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/ConstraintSetRemapping.cxxtest.hh
/// @brief  test suite for remapping the constraintSet upon conformation length changes
/// @author Florian Richter, Jan 09, floric@u.washington.edu

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <core/conformation/Residue.hh>


#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/scoring/constraints/BackboneStubConstraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

#include <numeric/conversions.hh>

//Auto Headers
#include <core/chemical/VariantType.hh>
#include <core/id/SequenceMapping.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/scoring/func/CircularHarmonicFunc.fwd.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.fwd.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.constraints.ConstraintSetRemapping.cxxtest");

using namespace core;

class ConstraintSetRemappingTests : public CxxTest::TestSuite
{

public:
	ConstraintSetRemappingTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	core::Real
	return_constraint_scores(  core::pose::Pose const & pose, utility::vector1< core::Size > const & positions ){

		using namespace::scoring;
		core::Real cst_score(0);

		for ( utility::vector1< core::Size >::const_iterator res_it = positions.begin(); res_it != positions.end(); ++res_it ) {
			EnergyMap scores = pose.energies().residue_total_energies( *res_it );

			cst_score += ( scores[ coordinate_constraint ] + scores[ atom_pair_constraint ] + scores[ angle_constraint ] + scores[ dihedral_constraint ] + scores[ res_type_constraint ] + scores[ backbone_stub_constraint ]);

		}

		return cst_score;
	}


	///////////////////////////////////////////////////////////////////////////////
	// ------------------------------------------ //
	/// @brief strategy: put in a couple of random constraints, remember the score, mess around with the pose length a bit,
	/// @brief then check if the score remains the same
	/// jan '11 added some tests to make sure that the constr
	void test_constraint_set_remapping()
	{
		using namespace core::scoring::constraints;
		using core::chemical::ResidueType;
		using core::chemical::AtomIndices;
		using core::conformation::Residue;
		using core::id::AtomID;

		pose::Pose pose( create_test_in_pdb_pose() );
		core::id::AtomID fixed_pt( pose.atom_tree().root()->atom_id() );
		//core::import_pose::pose_from_file( pose, "core/scoring/constraints/test_in.pdb" , core::import_pose::PDB_file);


		scoring::ScoreFunctionOP scorefxn( new scoring::ScoreFunction );
		scorefxn->reset();
		scorefxn->set_weight( scoring::fa_atr, 0.80 );
		scorefxn->set_weight( scoring::fa_rep, 0.44 );
		scorefxn->set_weight( scoring::fa_sol, 0.65 );
		scorefxn->set_weight( scoring::coordinate_constraint, 1.0 );
		scorefxn->set_weight( scoring::atom_pair_constraint, 1.0 );
		scorefxn->set_weight( scoring::angle_constraint, 1.0 );
		scorefxn->set_weight( scoring::dihedral_constraint, 1.0 );
		scorefxn->set_weight( scoring::res_type_constraint, 1.0 );
		scorefxn->set_weight( scoring::backbone_stub_constraint, 1.0 );

		core::Real const stddev_radians = numeric::conversions::radians( 5.0 );
		core::scoring::func::FuncOP some_func( new core::scoring::func::CircularHarmonicFunc( 1.05, stddev_radians ) );
		core::scoring::func::FuncOP some_other_func( new core::scoring::func::HarmonicFunc( 0.5, 1.0 ) );

		utility::vector1< core::Size > cst_positions;
		cst_positions.push_back( 3 );
		cst_positions.push_back( 12 );
		cst_positions.push_back( 15 );
		cst_positions.push_back( 17 );
		cst_positions.push_back( 19 );
		cst_positions.push_back( 20 );
		cst_positions.push_back( 30 );
		cst_positions.push_back( 33 );
		cst_positions.push_back( 35 );
		cst_positions.push_back( 44 );
		cst_positions.push_back( 45 );
		cst_positions.push_back( 48 );
		cst_positions.push_back( 57 );
		cst_positions.push_back( 65 );
		cst_positions.push_back( 69 );
		cst_positions.push_back( 70 );

		//some random constraints
		pose.add_constraint( scoring::constraints::ConstraintCOP( new AtomPairConstraint( AtomID(1, 12), AtomID(3, 57), some_func ) ) );
		pose.add_constraint( scoring::constraints::ConstraintCOP( new AtomPairConstraint( AtomID(6, 20), AtomID(1, 17), some_other_func ) ) );
		pose.add_constraint( scoring::constraints::ConstraintCOP( new AtomPairConstraint( AtomID(7, 15), AtomID(3, 16), some_other_func ) ) );
		pose.add_constraint( scoring::constraints::ConstraintCOP( new AngleConstraint( AtomID(6, 19), AtomID(1, 3), AtomID(10, 3 ), some_func ) ) );
		pose.add_constraint( scoring::constraints::ConstraintCOP( new AngleConstraint( AtomID(4, 69), AtomID(1, 57), AtomID(3, 30 ), some_func ) ) );
		pose.add_constraint( scoring::constraints::ConstraintCOP( new DihedralConstraint( AtomID(4, 65), AtomID(1, 65), AtomID(3, 65 ), AtomID(2, 65 ), some_func ) ) );
		pose.add_constraint( scoring::constraints::ConstraintCOP( new DihedralConstraint( AtomID(8, 44), AtomID(1, 44), AtomID(2, 48 ), AtomID(1, 48 ), some_func ) ) );
		pose.add_constraint( scoring::constraints::ConstraintCOP( new ResidueTypeConstraint( pose, 33, 1.0 ) ) );

		utility::vector1< ConstraintCOP > multi_csts;
		utility::vector1< ConstraintCOP > ambig_csts;

		multi_csts.push_back( ConstraintCOP( new AngleConstraint( AtomID(6, 15), AtomID(1, 15), AtomID(3, 3 ), some_func ) ) );
		multi_csts.push_back( ConstraintCOP( new AtomPairConstraint( AtomID(2, 15), AtomID(2, 3),  some_func ) ) );
		multi_csts.push_back( ConstraintCOP( new DihedralConstraint( AtomID(2, 15), AtomID(7, 15), AtomID(1, 65), AtomID(3, 65 ),  some_func ) ) );

		ambig_csts.push_back( ConstraintCOP( new AngleConstraint( AtomID(6, 33), AtomID(1, 33), AtomID(3, 30 ), some_func ) ) );
		ambig_csts.push_back( ConstraintCOP( new AtomPairConstraint( AtomID(2, 35), AtomID(2, 13),  some_func ) ) );
		ambig_csts.push_back( ConstraintCOP( new DihedralConstraint( AtomID(2, 35), AtomID(7, 35), AtomID(1, 45), AtomID(3, 65 ),  some_func ) ) );

		ConstraintCOP posemulticst = pose.add_constraint( scoring::constraints::ConstraintCOP( new MultiConstraint( multi_csts ) ) );
		ConstraintCOP poseambigcst = pose.add_constraint( scoring::constraints::ConstraintCOP( new AmbiguousConstraint( ambig_csts ) ) );

		//first test whether constraint == operators work
		utility::vector1< ConstraintCOP > testeq_csts;
		testeq_csts.push_back( ConstraintCOP( new AngleConstraint(AtomID(6, 15), AtomID(1, 15), AtomID(3, 3 ), some_func ) ) );
		testeq_csts.push_back( ConstraintCOP( new AtomPairConstraint( AtomID(2, 15), AtomID(2, 3),  some_func ) ) );
		testeq_csts.push_back( ConstraintCOP( new DihedralConstraint( AtomID(2, 15), AtomID(7, 15), AtomID(1, 65), AtomID(3, 65 ),  some_func ) ) );
		MultiConstraint testeq_mult( testeq_csts );
		AmbiguousConstraint testeq_ambig1( testeq_csts );

		utility::vector1< ConstraintCOP > testeq_ambig_csts;
		testeq_ambig_csts.push_back( ConstraintCOP( new AngleConstraint( AtomID(6, 33), AtomID(1, 33), AtomID(3, 30 ), some_func ) ) );
		testeq_ambig_csts.push_back( ConstraintCOP( new AtomPairConstraint( AtomID(2, 35), AtomID(2, 13),  some_func ) ) );
		testeq_ambig_csts.push_back( ConstraintCOP( new DihedralConstraint( AtomID(2, 35), AtomID(7, 35), AtomID(1, 45), AtomID(3, 65 ),  some_func ) ) );
		AmbiguousConstraint testeq_ambig2( testeq_ambig_csts );

		TS_ASSERT( *(testeq_csts[1]) == *(multi_csts[1]) );
		TS_ASSERT( *(testeq_csts[2]) == *(multi_csts[2]) );
		TS_ASSERT( *(testeq_csts[3]) == *(multi_csts[3]) );
		TS_ASSERT( *(testeq_csts[1]) != *(multi_csts[2]) );
		TS_ASSERT( *(testeq_csts[3]) != *(multi_csts[2]) );
		TS_ASSERT( testeq_mult == *posemulticst );
		TS_ASSERT( testeq_ambig1 != *posemulticst );
		TS_ASSERT( testeq_ambig2 == *poseambigcst );

		//equality operator tests done

		//score the pose and remeber the constraint score
		(*scorefxn)( pose );
		core::Real orig_score = return_constraint_scores( pose, cst_positions );

		TR << "orig_cst_score is  " << orig_score << std::endl;

		//now we go to town with the pose
		core::conformation::Residue dummy_ala( pose.residue( 3 ) );
		pose.append_polymer_residue_after_seqpos( dummy_ala, 19, true );
		for ( utility::vector1< core::Size >::iterator res_it = cst_positions.begin(); res_it != cst_positions.end(); ++res_it ) {
			if ( *res_it > 19 ) ( *res_it)++;
		}
		(*scorefxn)( pose );

		core::Real test1_score = return_constraint_scores( pose, cst_positions );
		TR << "test1_score (after appending at res19) is  " << test1_score << std::endl;
		TS_ASSERT_DELTA( orig_score, test1_score, 0.00001 );

		pose.delete_polymer_residue( 28 );
		pose.delete_polymer_residue( 28 );
		pose.delete_polymer_residue( 28 );
		for ( utility::vector1< core::Size >::iterator res_it = cst_positions.begin(); res_it != cst_positions.end(); ++res_it ) {
			if ( *res_it > 35 ) ( *res_it) -= 3;
		}
		(*scorefxn)( pose );

		core::Real test2_score = return_constraint_scores( pose, cst_positions );
		TR << "test2_score (after 3 deletions at pos 28) is  " << test2_score << std::endl;
		TS_ASSERT_DELTA( orig_score, test2_score, 0.00001 );

		pose.prepend_polymer_residue_before_seqpos( dummy_ala, 28, true );
		for ( utility::vector1< core::Size >::iterator res_it = cst_positions.begin(); res_it != cst_positions.end(); ++res_it ) {
			if ( *res_it > 28 ) ( *res_it)++;
		}
		(*scorefxn)( pose );

		core::Real test3_score = return_constraint_scores( pose, cst_positions );
		TR << "test3_score (after prepending before pos 28) is  " << test2_score << std::endl;
		TS_ASSERT_DELTA( orig_score, test3_score, 0.00001 );

		///now do some tests regarding adding of loop modelling residue types
		ConstraintCOP bstubcst =  pose.add_constraint( scoring::constraints::ConstraintCOP( new BackboneStubConstraint( pose, 15, fixed_pt, pose.residue(19), -80.0, 0.2) ) );
		(*scorefxn)( pose );
		core::Real test4_score = return_constraint_scores( pose, cst_positions );
		TR << "test4_score (after adding bstubcst at pos 15) is  " << test4_score << std::endl;
		pose::Pose tmppose = pose;
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, 15 );
		(*scorefxn)( pose );
		core::Real test5_score = return_constraint_scores( pose, cst_positions );
		TR << "test5_score (after adding cutpoint at pos 15) is  " << test5_score << std::endl;

		pose.constraint_set( tmppose.constraint_set()->remapped_clone( tmppose, pose ) );
		(*scorefxn)( pose );
		core::Real test6_score = return_constraint_scores( pose, cst_positions );
		TR << "test6_score (after adding cutpoint and remapping cstset at pos 15) is  " << test6_score << std::endl;
		TS_ASSERT_DELTA( test4_score, test6_score, 0.00001 );

		pose::Pose cutpose = pose;
		core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_LOWER, 15 );
		(*scorefxn)( pose );
		core::Real test7_score = return_constraint_scores( pose, cst_positions );
		TR << "test7_score (after removing cutpoint at pos 15) is  " << test7_score << std::endl;

		pose.constraint_set( cutpose.constraint_set()->remapped_clone( cutpose, pose ) );
		(*scorefxn)( pose );
		core::Real test8_score = return_constraint_scores( pose, cst_positions );
		TR << "test8_score (after removing cutpoint and remapping cstset at pos 15) is  " << test8_score << std::endl;
		TS_ASSERT_DELTA( test4_score, test8_score, 0.00001 );

		TS_ASSERT( !pose.remove_constraint( bstubcst ) );
		TS_ASSERT( pose.remove_constraint( bstubcst, true ) );

		(*scorefxn)( pose );
		core::Real test9_score = return_constraint_scores( pose, cst_positions );
		TR << "test9_score (after removing bstubcst at pos 15) is  " << test9_score << std::endl;
		TS_ASSERT_DELTA( test3_score, test9_score, 0.00001 );

	}


};
