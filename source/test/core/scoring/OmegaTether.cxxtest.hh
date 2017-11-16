// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
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

// Protocols headers -- to make it easier to build poses.
#include <protocols/cyclic_peptide/PeptideStubMover.hh>

//Minimizer
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/MoveMap.hh>

//Numeric headers
#include <numeric/angle.functions.hh>

using namespace std;

using core::Size;
using core::Real;
using core::pose::Pose;
using core::chemical::AA;


static basic::Tracer TR("core.scoring.OmegaTether.cxxtest");

class OmegaTetherTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
		pose_ = core::pose::PoseOP( new core::pose::Pose );
	}

	void tearDown() {
	}

	/// @brief Tests the omega scorefunction with alpha-amino acids.
	void test_alpha_aa_omega() {
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, "");
		builder.add_residue("Append", "ALA", 2, false, "N", 0, 1, "C");
		builder.add_residue("Append", "GLY:CtermProteinFull", 3, false, "N", 0, 2, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::omega, 0.5 );

		pose.set_omega(1, 180.0);

		utility::vector1< core::Real > energies_list;

		for ( core::Size j(0); j<360; j+=5 ) {
			pose.set_omega(2, static_cast<core::Real>(j));
			pose.update_residue_neighbors();

			energies_list.push_back( (*sfxn)(pose) );
			TR << j << "\t" << pose.omega(2) << "\t" << energies_list[energies_list.size()] << std::endl;

			//Delete the following: for debugging only
			//char outfile[256];
			//sprintf( outfile, "alpha_aa_omega_test_%03lu.pdb", j);
			//pose.dump_pdb(std::string(outfile));
		}

		core::Real last_energy(energies_list[energies_list.size()]);
		for ( core::Size i(1), imax(energies_list.size()); i<=imax; ++i ) {
			TS_ASSERT( std::abs( energies_list[i] - last_energy ) > 0.00001 ); //All energies different.
			last_energy = energies_list[i];
		}
	}

	/// @brief Tests the omega scorefunction with oligoureas.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_oligourea_omega() {
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, "");
		builder.add_residue("Append", "OU3_ALA", 2, false, "N", 0, 1, "C");
		builder.add_residue("Append", "GLY:CtermProteinFull", 3, false, "N", 0, 2, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::omega, 0.5 );

		pose.set_omega(1, 180.0);

		utility::vector1< core::Real > energies_list;

		for ( core::Size i(0); i<360; i+=5 ) {
			pose.set_mu( 2, static_cast<core::Real>(i)  );
			for ( core::Size j(0); j<360; j+=5 ) {
				pose.set_omega(2, static_cast<core::Real>(j));
				pose.update_residue_neighbors();

				energies_list.push_back( (*sfxn)(pose) );
				TR << i << "\t" << j << "\t" << pose.mu(2) << "\t" << pose.omega(2) << "\t" << energies_list[energies_list.size()] << std::endl;

				//Delete the following: for debugging only
				//char outfile[256];
				//sprintf( outfile, "oligourea_omega_test_%03lu_%03lu.pdb", i, j);
				//pose.dump_pdb(std::string(outfile));
			}
		}
		core::Real last_energy(energies_list[energies_list.size()]);
		for ( core::Size i(1), imax(energies_list.size()); i<=imax; ++i ) {
			if ( i!=73 ) { TS_ASSERT( std::abs( energies_list[i] - last_energy ) > 0.00001 ); } //All energies different.
			else { TS_ASSERT_DELTA( energies_list[i], last_energy, 0.00001 ); } //These two are actually the same.
			last_energy = energies_list[i];
		}
	}

	/// @brief Tests the omega scorefunction minimization with oligoureas.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_oligourea_omega_min() {
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, "");
		builder.add_residue("Append", "OU3_ALA", 2, false, "N", 0, 1, "C");
		builder.add_residue("Append", "GLY:CtermProteinFull", 3, false, "N", 0, 2, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::omega, 0.5 );

		pose.set_omega(1, 180.0);

		pose.set_omega(2, 173.0);
		pose.set_mu(2, 168.0);

		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		for ( core::Size i=1; i<=3; ++i )  {
			mm->set_bb(i, true);
			mm->set_chi(i, false);
		}
		core::optimization::AtomTreeMinimizer minimizer;
		core::optimization::MinimizerOptionsOP min_options( new core::optimization::MinimizerOptions( "linmin_iterated", 0.0000001, true, false, false ) );
		minimizer.run( pose, *mm, *sfxn, *min_options );

		TS_ASSERT_DELTA( numeric::nonnegative_principal_angle_degrees( pose.omega(2) ), 180.0, 0.5 );
		TS_ASSERT_DELTA( numeric::nonnegative_principal_angle_degrees( pose.mu(2) ), 180.0, 0.5 );

		pose.set_omega(2, -173.0);
		pose.set_mu(2, -168.0);

		minimizer.run(pose, *mm, *sfxn, *min_options);

		TS_ASSERT_DELTA( numeric::nonnegative_principal_angle_degrees( pose.omega(2) ), 180.0, 0.5 );
		TS_ASSERT_DELTA( numeric::nonnegative_principal_angle_degrees( pose.mu(2) ), 180.0, 0.5 );
	}

	/// @brief Tests minimization of beta-amino acids with the omega scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_beta_aa_omega_min() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::omega, 0.5 );

		//Set up the pose
		pose_->clear();
		pose_ = core::import_pose::pose_from_file("core/scoring/betapose.pdb", core::import_pose::PDB_file);
		pose_->set_phi(3, 140.0);
		pose_->set_theta(3, /*140.0*/134.448); // new optimum for beta_nov15
		pose_->set_psi(3, /*140.0*/103.883);// new optimum for beta_nov15
		pose_->set_omega(3, 140.0);
		pose_->dump_pdb("core/scoring/betapose.premin.pdb");

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

		pose_->dump_pdb("core/scoring/betapose.min.pdb");

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
