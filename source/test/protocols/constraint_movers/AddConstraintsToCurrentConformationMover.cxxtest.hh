// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/constraint_movers/AddConstraintsToCurrentConformationMover.cxxtest.hh
/// @brief  test for AddConstraintsToCurrentConformationMover mover
/// @author Steven Lewis smlewi@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/conformation/Residue.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/constraint_movers/AddConstraintsToCurrentConformationMover.hh>

#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/Func.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <memory>

#include <core/scoring/ScoreFunction.hh> // AUTO IWYU For ScoreFunction
#include <core/scoring/ScoreFunctionFactory.hh> // AUTO IWYU For ScoreFunctionFactory

static basic::Tracer TR("protocols.constraint_movers.AddConstraintsToCurrentConformationMover.cxxtest.hh");

// --------------- Test Class --------------- //

class AddConstraintsToCurrentConformationMoverTests : public CxxTest::TestSuite {

private:
	core::pose::PoseOP test_dimer_pose_;
	core::scoring::ScoreFunctionOP scorefxn_;
public:

	void setUp() {
		core_init();

		test_dimer_pose_ = create_2res_1ten_2res_trp_cage_poseop(); //dimer structure
	}

	void tearDown() {
	}

	// //this was used for setting the tests up
	// void cst_shower() {
	//
	//  core::scoring::constraints::ConstraintSetCOP cstset( test_dimer_pose_->constraint_set());
	//
	//  //core::scoring::constraints::ConstraintCOPs
	//  auto csts(cstset->get_all_constraints());
	//  TR << csts.size() << std::endl;
	//  for (auto eachcst : csts){
	//   TR << eachcst->type() << std::endl; //show(TR);
	//   eachcst->get_func().show_definition(TR);
	//   TR << std::endl;
	//  }
	// }

	//each of these tests tests a similar thing
	void multi_test( core::Size const num_cst, std::string const & cst_type, std::string const & func_type) {

		core::scoring::constraints::ConstraintSetCOP cstset( test_dimer_pose_->constraint_set());

		TS_ASSERT(!cstset->is_empty());
		TS_ASSERT(cstset->has_constraints());
		TS_ASSERT(cstset->has_residue_pair_constraints());

		//core::scoring::constraints::ConstraintCOPs
		auto csts(cstset->get_all_constraints());
		TS_ASSERT_EQUALS(csts.size(), num_cst);  //how many constraints?
		for ( auto eachcst : csts ) {
			TS_ASSERT_EQUALS(cst_type, eachcst->type()); //what type of constraint?
			std::ostringstream func_def;
			eachcst->get_func().show_definition(func_def);
			//split definition into vector of strings; first element is func type
			//what type of func?
			TS_ASSERT_EQUALS(func_type, utility::split(func_def.str())[1]);
		}
		return;
	}


	//SML Dec 19 2016 - this is obviously not a complete test but it's better than zero

	void test_atompair_bounded() {
		protocols::constraint_movers::AddConstraintsToCurrentConformationMover actccm;
		actccm.use_distance_cst() = true;
		actccm.use_bounded_func() = true;
		actccm.bound_width() = 8.0;
		actccm.max_distance() = 10000; //just something large, I have no idea how far apart the test pose residues are
		actccm.min_seq_sep() = 1;

		actccm.apply(*test_dimer_pose_);
		TR << "TEST 1" << std::endl;
		//test_dimer_pose_->constraint_set()->show(TR);
		//TR << std::endl;
		//cst_shower();
		multi_test(6, "AtomPair", "BOUNDED");

	}

	void test_atompair_scalar_weighted() {
		protocols::constraint_movers::AddConstraintsToCurrentConformationMover actccm;
		actccm.use_distance_cst() = true;
		actccm.use_bounded_func() = false;
		actccm.coord_dev() = 1.0;
		actccm.max_distance() = 10000; //just something large, I have no idea how far apart the test pose residues are
		actccm.min_seq_sep() = 1;

		actccm.apply(*test_dimer_pose_);
		TR << "TEST 2" << std::endl;
		//test_dimer_pose_->constraint_set()->show(TR);
		//TR << std::endl;
		//cst_shower();
		multi_test(6, "AtomPair", "SCALARWEIGHTEDFUNC");

	}


	void test_atompair_harmonic() {
		protocols::constraint_movers::AddConstraintsToCurrentConformationMover actccm;
		actccm.use_distance_cst() = true;
		actccm.use_bounded_func() = false;
		actccm.use_harmonic_func() = true;
		actccm.coord_dev() = 1.0;
		actccm.min_seq_sep() = 1;
		actccm.max_distance() = 10000; // just something large, we have no idea how far apart the test pose residues are

		actccm.apply(*test_dimer_pose_);
		multi_test(6, "AtomPair", "HARMONIC");

	}

	void test_coord_bounded() {
		protocols::constraint_movers::AddConstraintsToCurrentConformationMover actccm;
		actccm.use_distance_cst() = false;
		actccm.use_bounded_func() = true;
		actccm.bound_width() = 8.0;

		actccm.apply(*test_dimer_pose_);
		//TR << "TEST 3" << std::endl;
		//test_dimer_pose_->constraint_set()->show(TR);
		//TR << std::endl;
		//cst_shower();
		multi_test(4, "CoordinateConstraint", "BOUNDED");

	}

	void test_coord_harmonic() {
		protocols::constraint_movers::AddConstraintsToCurrentConformationMover actccm;
		actccm.use_distance_cst() = false;
		actccm.use_bounded_func() = false;
		actccm.coord_dev() = 1.0;

		actccm.apply(*test_dimer_pose_);
		//TR << "TEST 4" << std::endl;
		//test_dimer_pose_->constraint_set()->show(TR);
		//TR << std::endl;
		//cst_shower();
		multi_test(4, "CoordinateConstraint", "HARMONIC");

	}

	void helper__test_coord_recovery(){

	}

	void test_sc_tip(){
		using namespace protocols::constraint_movers;
		AddConstraintsToCurrentConformationMover actccm;
		actccm.use_distance_cst() = false;
		actccm.set_atom_selector( AddConstraintsToCurrentConformationMover::AtomSelector::SC_TIP_ONLY );

		using namespace core::scoring;
		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "none.wts" );
		sfxn->set_weight( coordinate_constraint, 100.0 );

		core::pose::Pose pose;
		core::pose::make_pose_from_sequence( pose, "E/E", "fa_standard" );
		TS_ASSERT_EQUALS( pose.num_chains(), 2 );

		core::Size const cst_resid = 2; //we'll be moving this one later

		// STORE STARTING XYZs FOR RES 2:
		utility::vector1< numeric::xyzVector< core::Real > > starting_xyzs;
		for ( core::Size a = 1; a <= pose.residue( cst_resid ).nheavyatoms(); ++a ) {
			starting_xyzs.push_back( pose.residue( cst_resid ).xyz( a ) );
		}

		// APPLY CSTS AND ASSERT LOW SCORE
		// (score should be low because nothing has moved yet
		actccm.apply( pose );
		TS_ASSERT_DELTA( sfxn->score( pose ), 0.0, 0.001 );

		// ANALYZE CSTS
		core::scoring::constraints::ConstraintCOPs const & csts =
			pose.constraint_set()->get_all_constraints();
		core::Size n_csts_for_resid = 0;
		for ( core::scoring::constraints::ConstraintCOP cst : csts ) {
			core::scoring::constraints::CoordinateConstraintCOP cc =
				std::dynamic_pointer_cast< core::scoring::constraints::CoordinateConstraint const >( cst );
			auto const atomid = cc->atom( 1 );
			if ( atomid.rsd() == cst_resid ) {
				++n_csts_for_resid;

				// ASSERT THAT THE CST XYZs MATCH OUR STARTING STATE
				auto const & vec = cc->xyz_target();
				TS_ASSERT_DELTA( starting_xyzs[ atomid.atomno() ].x(), vec.x(), 0.0001 );
				TS_ASSERT_DELTA( starting_xyzs[ atomid.atomno() ].y(), vec.y(), 0.0001 );
				TS_ASSERT_DELTA( starting_xyzs[ atomid.atomno() ].z(), vec.z(), 0.0001 );
			}
		}

		// ASSERT THAT WE HAVE THE NUMBER OF CSTS WE EXPECT
		constexpr core::Size atoms_per_E = 4; //only the 4 atoms at the tip of E should be included
		constexpr core::Size nE = 2; //We have 2 E residues in this pose
		TS_ASSERT_EQUALS( n_csts_for_resid, atoms_per_E );
		TS_ASSERT_EQUALS( csts.size(), atoms_per_E*nE );

		//MINIMIZE AND ASSERT THAT NOTHING MOVES
		core::select::movemap::MoveMapFactory mmf;
		mmf.all_bb( false );
		mmf.all_chi( false );
		mmf.all_jumps( true );

		using namespace protocols::minimization_packing;
		MinMover min;
		min.score_function( sfxn );
		min.movemap_factory( mmf.clone() );
		min.tolerance( 0.000001 );
		min.max_iter( 999999999 );
		min.min_type( "dfpmin" );//No hill-climbing according to VKM

		min.apply( pose );

		TS_ASSERT_DELTA( sfxn->score( pose ), 0.0, 0.001 );

		// MOVER THE POSE JUST A LITTLE
		// Chi 3 += 15 degrees
		pose.set_chi( 3, cst_resid, pose.chi( 3, cst_resid ) + 15.0 );
		//Both of these two atoms should move
		core::Size const atm1 = starting_xyzs.size() - 1;
		core::Size const atm2 = starting_xyzs.size();

		// ASSERT THAT THE TWO ATOMS MOVE MORE THAN 0.1 ANGSTROMS
		TS_ASSERT( starting_xyzs[ atm1 ].distance( pose.residue( cst_resid ).xyz( atm1 ) ) > 0.1 );
		TS_ASSERT( starting_xyzs[ atm2 ].distance( pose.residue( cst_resid ).xyz( atm2 ) ) > 0.1 );
		TS_ASSERT( sfxn->score( pose ) > 0.0 );
		core::Real const bad_score = sfxn->score( pose );
		TR << "Score after chi move: " << bad_score << std::endl;

		// REVERT TO ORIGINAL STATE
		min.apply( pose );
		TR << "Score after min: " << sfxn->score( pose ) << std::endl;

		// ASSERT REVERSION
		// Score might not go down to 0 but should be close and better than before
		TS_ASSERT( sfxn->score( pose ) < bad_score );
		TS_ASSERT( starting_xyzs[ atm1 ].distance( pose.residue( cst_resid ).xyz( atm1 ) ) < 0.1 );//THIS FAILS
		TS_ASSERT( starting_xyzs[ atm2 ].distance( pose.residue( cst_resid ).xyz( atm2 ) ) < 0.1 );//THIS FAILS
		TS_ASSERT_DELTA( sfxn->score( pose ), 0.0, 0.01 );
	}

};
