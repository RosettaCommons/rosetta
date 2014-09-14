// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/generalized_kinematic_closure/GeneralizedKIC.cxxtest.hh
/// @brief  Unit tests for the GeneralizedKIC mover.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers:
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pdb1ubq.hh>

// GeneralizedKIC headers:
#include <protocols/generalized_kinematic_closure/GeneralizedKIC.hh>
#include <protocols/generalized_kinematic_closure/util.hh>

// Other Rosetta libraries:
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/simple_moves/MutateResidue.hh>


// --------------- Test Class --------------- //

class GeneralizedKIC_Tests : public CxxTest::TestSuite {

private:
	core::pose::PoseOP testpose_;
	core::scoring::ScoreFunctionOP scorefxn_;

public:

	void setUp() {
		core_init_with_additional_options( "-extra_res_fa protocols/generalized_kinematic_closure/DALA.params" );

		scorefxn_ = core::scoring::get_score_function();
		testpose_ = pdb1ubq5to13_poseop();

	}

	void tearDown() {
	}

	/// @brief Test GenKIC as applied to an L-alpha backbone loop.
	/// @details This uses residues 5-13 from the 1ubq structure, and
	/// samples conformations for resiudes 6-12.  It then double-checks
  /// that (a) at least one solution was found, and (b) that all
  /// peptide bond lengths are the same as in the original structure.
	void test_GeneralizedKIC_L_alpha_backbone()
	{
		using namespace protocols::generalized_kinematic_closure;

		GeneralizedKICOP genkic( new GeneralizedKIC ); //Create the mover.

		//Add the loop residues:
		for( core::Size i=2; i<=8; ++i ) genkic->add_loop_residue(i);

		//Set the pivots:
		genkic->set_pivot_atoms( 2, "CA", 4, "CA", 8, "CA" );

		//Add a randomizing-by-rama perturber:
		genkic->add_perturber( "randomize_alpha_backbone_by_rama" );
		for( core::Size i=2; i<=8; ++i ) genkic->add_residue_to_perturber_residue_list( i );

		//Add a loop bump check:
		genkic->add_filter( "loop_bump_check" );

		//Set options:
		genkic->set_closure_attempts(25000); //Try a maximum of 25000 times.
		genkic->set_n_closure_attempts_is_a_maximum(true); //Stop when a solution is found.

		//Add a lowest_energy selector:
		genkic->set_selector_type("lowest_energy_selector");
		genkic->set_selector_scorefunction(scorefxn_);


		core::pose::PoseOP testpose2 = testpose_->clone();

		genkic->apply(*testpose2);
		//testpose_->dump_pdb("temp1.pdb"); //DELETE ME
		//testpose2->dump_pdb("temp2.pdb"); //DELETE ME

		//CHECK A: Was the closure successful?
		TS_ASSERT(genkic->last_run_successful());

		//CHECK B: Are all the bond lengths as before?
		for(core::Size i=1; i<9; ++i) { //Loop through all residues
			core::Real peplength1 = testpose_->residue(i).xyz("C").distance( testpose_->residue(i+1).xyz("N") );
			core::Real peplength2 = testpose2->residue(i).xyz("C").distance( testpose2->residue(i+1).xyz("N") );
			TS_ASSERT_DELTA( peplength1, peplength2, 1.0e-3 );
		}

		return;
	}

	/// @brief Test GenKIC as applied to an alpha backbone loop containing a D-residue.
	/// @details This uses residues 5-13 from the 1ubq structure with G10 mutated to
	/// DALA, and samples conformations for resiudes 6-12.  It then double-checks
  /// that (a) at least one solution was found, and (b) that all
  /// peptide bond lengths are the same as in the original structure.
	void test_GeneralizedKIC_mixed_DL_alpha_backbone()
	{
		using namespace protocols::generalized_kinematic_closure;

		GeneralizedKICOP genkic( new GeneralizedKIC ); //Create the mover.


		//Copy the pose, and mutate a residue to DALA:
		core::pose::PoseOP testpose2 = testpose_->clone();
		protocols::simple_moves::MutateResidue mutres(6, "DALA");
		mutres.apply(*testpose2);
		
		//Add the loop residues:
		for( core::Size i=2; i<=8; ++i ) genkic->add_loop_residue(i);

		//Set the pivots:
		genkic->set_pivot_atoms( 2, "CA", 4, "CA", 8, "CA" );

		//Add a randomizing-by-rama perturber:
		genkic->add_perturber( "randomize_alpha_backbone_by_rama" );
		for( core::Size i=2; i<=8; ++i ) genkic->add_residue_to_perturber_residue_list( i );

		//Add a loop bump check:
		genkic->add_filter( "loop_bump_check" );

		//Set options:
		genkic->set_closure_attempts(25000); //Try a maximum of 25000 times.
		genkic->set_n_closure_attempts_is_a_maximum(true); //Stop when a solution is found.

		//Add a lowest_energy selector:
		genkic->set_selector_type("lowest_energy_selector");
		genkic->set_selector_scorefunction(scorefxn_);

		genkic->apply(*testpose2);
		//testpose_->dump_pdb("temp3.pdb"); //DELETE ME
		//testpose2->dump_pdb("temp4.pdb"); //DELETE ME

		//CHECK A: Was the closure successful?
		TS_ASSERT(genkic->last_run_successful());

		//CHECK B: Are all the bond lengths as before?
		for(core::Size i=1; i<9; ++i) { //Loop through all residues
			core::Real peplength1 = testpose_->residue(i).xyz("C").distance( testpose_->residue(i+1).xyz("N") );
			core::Real peplength2 = testpose2->residue(i).xyz("C").distance( testpose2->residue(i+1).xyz("N") );
			TS_ASSERT_DELTA( peplength1, peplength2, 1.0e-3 );
		}

		return;
	}


}; //class GeneralizedKIC_Tests
