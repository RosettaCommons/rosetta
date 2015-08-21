// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file test/core/scoring/OmegaTether.cxxtest.hh
/// @brief Unit tests for the omega score term.
/// @detials Gaah!  The omega scoreterm is not covered by unit tests.  I've only added one for beta-amino acids.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

// Unit headers
#include <core/scoring/OmegaTether.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunction.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/AA.hh>

//Minimizer
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/MoveMap.hh>

using namespace std;

using core::Size;
using core::Real;
using core::pose::Pose;
using core::chemical::AA;


static thread_local basic::Tracer TR("core.scoring.OmegaTether.cxxtest");

class OmegaTetherTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init_with_additional_options( "-extra_res_fa beta-peptide/B3A.params" );
		pose_ = core::pose::PoseOP( new core::pose::Pose );
	}

	void tearDown() {
	}

	/// @brief Tests minimization of beta-amino acids with the omega scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_beta_aa_omega_min() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::omega, 0.5 );

		//Set up the pose
		pose_->clear();
		pose_ = core::import_pose::pose_from_pdb("core/scoring/betapose.pdb");
		pose_->set_phi(3, 140.0);
		pose_->set_theta(3, 140.0);
		pose_->set_psi(3, 140.0);
		pose_->set_omega(3, 140.0);

		core::Real const phi_start( pose_->phi(3));
		core::Real const theta_start( pose_->theta(3));
		core::Real const psi_start( pose_->psi(3));
		core::Real const omega_start( pose_->omega(3));

		//Set up minimizer
		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		for ( core::Size i=1; i<=5; ++i ) {
			mm->set_bb(i,true);
			mm->set_chi(i,false);
		}

		core::optimization::AtomTreeMinimizer minimizer;
		core::optimization::MinimizerOptionsOP min_options( new core::optimization::MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.0000001, true, false, false ) );

		minimizer.run( *pose_, *mm, *scorefxn, *min_options );

		core::Real const phi_end( pose_->phi(3) );
		core::Real const theta_end( pose_->theta(3) );
		core::Real const psi_end( pose_->psi(3) );
		core::Real const omega_end( pose_->omega(3) );

		if ( TR.visible() ) {
			TR << "ANGLE\tSTART\tEND" << std::endl;
			TR << "Phi:\t" << phi_start << "\t" << phi_end << std::endl;
			TR << "Theta:\t" << theta_start << "\t" << theta_end << std::endl;
			TR << "Psi:\t" << psi_start << "\t" << psi_end << std::endl;
			TR << "Omega:\t" << omega_start << "\t" << omega_end << std::endl;
		}

		//None of these should change much, if at all:
		TS_ASSERT_DELTA( phi_start, phi_end, 0.1 );
		TS_ASSERT_DELTA( theta_start, theta_end, 0.1 );
		TS_ASSERT_DELTA( psi_start, psi_end, 0.1 );

		//This should change more (starts at about 140):
		TS_ASSERT_DELTA(omega_end, 180.0, 0.5);

		pose_->clear();
	}

private:
	core::pose::PoseOP pose_;

};
