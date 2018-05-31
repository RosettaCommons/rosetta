// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file test/core/scoring/SymmGlyMin.cxxtest.hh
/// @brief Unit tests for glycine minimization with the --symmetric_gly_tables option.
/// @detials Left- and right-handed conformations of glycine should minimize identically with this option.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

// Unit headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/variant_util.hh>
#include <core/chemical/AA.hh>

//Minimizer
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/MoveMap.hh>

using namespace std;

using core::Size;
using core::Real;
using core::pose::Pose;
using core::chemical::AA;


static basic::Tracer TR("core.scoring.SymmGlyMinTests.cxxtest");

class SymmGlyMinTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init_with_additional_options( "-symmetric_gly_tables true -run:constant_seed -nodelay -run:jran 35153" );
	}

	void tearDown() {
	}

	/// @brief Run the minimizer on the pose.
	///
	void do_minimization( core::pose::Pose &pose, core::scoring::ScoreFunctionOP sfxn, bool const cartesian ) {
		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		mm->set_bb(true);
		mm->set_chi(true);
		if ( cartesian ) {
			core::optimization::CartesianMinimizer minimizer;
			core::optimization::MinimizerOptionsOP min_options( new core::optimization::MinimizerOptions( "linmin", 10.0, true, false, false ) );
			minimizer.run( pose, *mm, *sfxn, *min_options );
		} else {
			core::optimization::AtomTreeMinimizer minimizer;
			core::optimization::MinimizerOptionsOP min_options( new core::optimization::MinimizerOptions( "linmin", 10.0, true, false, false ) );
			minimizer.run( pose, *mm, *sfxn, *min_options );
		}
	}

	/// @brief Are two angles within a threshhold of one another?
	///
	bool within_thresh( core::Real const &val1, core::Real const &val2, core::Real const &thresh ) {
		if ( std::abs( val1 - val2 ) <= thresh ) return true;
		core::Real val1prime = val1;
		core::Real val2prime = val2;
		if ( val1<val2 ) {
			val1prime += 360.0;
		} else {
			val2prime += 360.0;
		}
		return ( std::abs( val1prime - val2prime ) <= thresh );
	}

	/// @brief Construct repeat sequences with poly-glycine, and confirm that mirror-image conformations
	/// minimize identically with a given scorefunction.
	void repeat_structure_test( core::scoring::ScoreFunctionOP sfxn, bool const cartesian ) {
		//int count = 0; //DELETE ME

		core::pose::PoseOP pose( new core::pose::Pose() );
		core::pose::make_pose_from_sequence(*pose, "GGGGGGGG", "fa_standard", cartesian);
		core::pose::PoseOP pose2( pose->clone() );
		for ( int iphi=-180; iphi<180; iphi+=60 ) {
			for ( int ipsi=-180; ipsi<=180; ipsi+=60 ) {
				core::pose::PoseOP pose( new core::pose::Pose() );
				core::pose::make_pose_from_sequence(*pose, "GGGGGGGG", "fa_standard", false);
				core::pose::PoseOP pose2( pose->clone() );

				//The following lines may be deleted:
				/*++count;
				char fname1[256];
				char fname2[256];
				char fname3[256];
				sprintf(fname1, "vsymmglytest_pose1_%04i.pdb", count);
				sprintf(fname2, "vsymmglytest_pose2_%04i.pdb", count);
				sprintf(fname3, "vsymmglytest_pose3_%04i.pdb", count);*/

				for ( core::Size ir=1; ir<=8; ++ir ) {
					pose->set_omega(ir, 180.0);
					pose->set_phi(ir, static_cast<core::Real>(iphi));
					pose->set_psi(ir, static_cast<core::Real>(ipsi));
					pose2->set_omega(ir, 180.0);
					pose2->set_phi(ir, -1.0*static_cast<core::Real>(iphi));
					pose2->set_psi(ir, -1.0*static_cast<core::Real>(ipsi));
				}

				//Make a third pose, too, and mirror this one:
				core::pose::PoseOP pose3( pose->clone() );
				for ( core::Size ir=1, irmax=pose->total_residue(); ir<=irmax; ++ir ) {
					core::conformation::ResidueOP ires( pose->residue( ir ).clone() );
					ires->set_mirrored_relative_to_type(true);
					for ( core::Size ia=1, iamax=pose->residue_type(ir).natoms(); ia<=iamax; ++ia ) {
						core::id::AtomID const curat( ia, ir );
						numeric::xyzVector< core::Real > v( pose->xyz(curat) );
						v.x() *= -1.0;
						ires->set_xyz( ia, v );
					}
					pose3->replace_residue( ir, *ires, false );
				}
				pose3->update_residue_neighbors();

				core::pose::PoseOP pose4a, pose4b;
				core::pose::PoseOP pose5a, pose5b;
				if ( cartesian ) { //Additional poses used for additional cartesian tests.
					pose4a = pose->clone();
					pose4b = pose3->clone();
					pose5a = pose->clone();
					pose5b = pose2->clone();
					core::pose::add_lower_terminus_type_to_pose_residue(*pose5a, 1);
					core::pose::add_lower_terminus_type_to_pose_residue(*pose5b, 1);
					core::pose::add_upper_terminus_type_to_pose_residue(*pose5a, pose5a->total_residue());
					core::pose::add_upper_terminus_type_to_pose_residue(*pose5b, pose5b->total_residue());
				}

				(*sfxn)(*pose);
				do_minimization(*pose, sfxn, cartesian);
				(*sfxn)(*pose2);
				do_minimization(*pose2, sfxn, cartesian);
				(*sfxn)(*pose3);
				do_minimization(*pose3, sfxn, cartesian);

				//Delete the following:
				/*pose->dump_pdb(std::string(fname1));
				pose2->dump_pdb(std::string(fname2));
				pose3->dump_pdb(std::string(fname3));*/

				bool const skip_pose3_tests( sfxn->get_weight( core::scoring::rama ) && ( iphi == 0 || ipsi == 0 ) ); //The phi=0 || psi=0 lines have a known Rama discontinuity.  Blargh.

				TS_ASSERT_DELTA(pose->energies().total_energy(), pose2->energies().total_energy(), std::max( std::abs( std::max(pose->energies().total_energy(), pose2->energies().total_energy())/1000.0 ), 1e-12 ) );
				if ( !skip_pose3_tests ) TS_ASSERT_DELTA(pose->energies().total_energy(), pose3->energies().total_energy(), std::max( std::abs( std::max(pose->energies().total_energy(), pose3->energies().total_energy())/1000.0 ), 1e-12 ) ); //Skip this check for phi=0, psi=0 and the rama score term.  There's a known silly discontinuity there.1
				TR << "E1\t" << pose->energies().total_energy() << "\tE2\t" << pose2->energies().total_energy() << "\tE3\t" << pose3->energies().total_energy() << std::endl;
				for ( core::Size ir=1, irmax=pose->size(); ir<=irmax; ++ir ) {
					TR << "iphi\t" << iphi << "\tipsi\t" << ipsi << std::endl;
					TR << "phi1\t" << pose->phi(ir) << "\tphi2\t" << pose2->phi(ir) << "\tphi3\t" << pose3->phi(ir) << std::endl;
					TR << "psi1\t" << pose->psi(ir) << "\tpsi2\t" << pose2->psi(ir) << "\tpsi3\t" << pose3->psi(ir) << std::endl;
					TR << "omega1\t" << pose->omega(ir) << "\tomega2\t" << pose2->omega(ir) << "\tomega3\t" << pose3->omega(ir) << std::endl;
					TS_ASSERT( within_thresh( pose->phi(ir), -1.0*pose2->phi(ir), 0.01 ) );
					if ( !skip_pose3_tests ) TS_ASSERT( within_thresh( pose->phi(ir), -1.0*pose3->phi(ir), 0.01 ) );
					TS_ASSERT( within_thresh( pose->psi(ir), -1.0*pose2->psi(ir), 0.01 ) );
					if ( !skip_pose3_tests ) TS_ASSERT( within_thresh( pose->psi(ir), -1.0*pose3->psi(ir), 0.01 ) );
					TS_ASSERT( within_thresh( pose->omega(ir), -1.0*pose2->omega(ir), 0.01 ) );
					if ( !skip_pose3_tests ) TS_ASSERT( within_thresh( pose->omega(ir), -1.0*pose3->omega(ir), 0.01 ) );
				}

				if ( cartesian ) { //For cartesian tests, do another test where we jitter the poses.
					//Jitter coordinates of pose4a, and mirror them in pose4b:
					for ( core::Size ir=1, irmax=pose4a->total_residue(); ir<=irmax; ++ir ) {
						for ( core::Size ia=1, iamax=pose4a->residue_type(ir).natoms(); ia<=iamax; ++ia ) {
							core::id::AtomID const curat( ia, ir );
							numeric::xyzVector< core::Real > curpos( pose4a->xyz( curat ) );
							curpos.x() += numeric::random::rg().gaussian() * 3.0;
							curpos.y() += numeric::random::rg().gaussian() * 3.0;
							curpos.z() += numeric::random::rg().gaussian() * 3.0;
							pose4a->set_xyz( curat, curpos );
							curpos.x() *= -1.0;
							pose4b->set_xyz( curat, curpos );
						}
					}
					(*sfxn)(*pose4a);
					(*sfxn)(*pose4b);
					TR << "E4a_pre\t" << pose4a->energies().total_energy() << "\tE4b_pre\t" << pose4b->energies().total_energy() << std::endl;
					if ( !skip_pose3_tests ) TS_ASSERT_DELTA(pose4a->energies().total_energy(), pose4b->energies().total_energy(), std::max( std::abs( std::max(pose4a->energies().total_energy(), pose4b->energies().total_energy())/1000.0 ), 1e-12 ) );
					do_minimization(*pose4a, sfxn, cartesian);
					do_minimization(*pose4b, sfxn, cartesian);
					TR << "E4a_post\t" << pose4a->energies().total_energy() << "\tE4b_post\t" << pose4b->energies().total_energy() << std::endl;
					if ( !skip_pose3_tests ) TS_ASSERT_DELTA(pose4a->energies().total_energy(), pose4b->energies().total_energy(), std::max( std::abs( std::max(pose4a->energies().total_energy(), pose4b->energies().total_energy())/1000.0 ), 1e-12 ) );
					for ( core::Size ir=1, irmax=pose4a->size(); ir<=irmax; ++ir ) {
						TR << "iphi\t" << iphi << "\tipsi\t" << ipsi << std::endl;
						TR << "phi4a\t" << pose4a->phi(ir) << "\tphi4b\t" << pose4b->phi(ir) << std::endl;
						TR << "psi4a\t" << pose4a->psi(ir) << "\tpsi4b\t" << pose4b->psi(ir) << std::endl;
						TR << "omega4a\t" << pose4a->omega(ir) << "\tomega4b\t" << pose4b->omega(ir) << std::endl;
						if ( !skip_pose3_tests ) {
							TS_ASSERT( within_thresh( pose4a->phi(ir), -1.0*pose4b->phi(ir), 0.01 ) );
							TS_ASSERT( within_thresh( pose4a->psi(ir), -1.0*pose4b->psi(ir), 0.01 ) );
							TS_ASSERT( within_thresh( pose4a->omega(ir), -1.0*pose4b->omega(ir), 0.01 ) );
						}
					}

					//Testing torsion potentials for omega:
					for ( core::Size ir=1, irmax=pose5a->total_residue(); ir<=irmax; ++ir ) {
						pose5a->set_omega(ir, static_cast<core::Real>(iphi));
						pose5b->set_omega(ir, -1.0*static_cast<core::Real>(iphi));
					}
					pose5a->update_residue_neighbors();
					pose5b->update_residue_neighbors();
					(*sfxn)(*pose5a);
					(*sfxn)(*pose5b);
					TR << "E5a_pre\t" << pose5a->energies().total_energy() << "\tE5b_pre\t" << pose5b->energies().total_energy() << std::endl;
					if ( !skip_pose3_tests ) {
						TS_ASSERT_DELTA(pose5a->energies().total_energy(), pose5b->energies().total_energy(), std::max( std::abs( std::max(pose5a->energies().total_energy(), pose5b->energies().total_energy())/1000.0 ), 1e-12 ) );
						if ( std::abs( pose5a->energies().total_energy() - pose5b->energies().total_energy() ) > std::max( std::abs( std::max(pose5a->energies().total_energy(), pose5b->energies().total_energy())/1000.0 ), 1e-12 ) ) {
							TR << "Failure on phi=" << iphi << " psi=" << ipsi << " omega=" << iphi << ".  Dumping..." << std::endl;
							char outstr[256];
							sprintf( outstr, "pre_failure_A_%03i_%03i.pdb", iphi, ipsi );
							pose5a->dump_pdb(std::string(outstr));
							sprintf( outstr, "pre_failure_B_%03i_%03i.pdb", iphi, ipsi );
							pose5b->dump_pdb(std::string(outstr));
						}
					}
					do_minimization(*pose5a, sfxn, cartesian);
					do_minimization(*pose5b, sfxn, cartesian);
					pose5a->update_residue_neighbors();
					pose5b->update_residue_neighbors();
					(*sfxn)(*pose5a);
					(*sfxn)(*pose5b);
					TR << "E5a_post\t" << pose5a->energies().total_energy() << "\tE5b_post\t" << pose5b->energies().total_energy() << std::endl;
					if ( !skip_pose3_tests ) {
						TS_ASSERT_DELTA(pose5a->energies().total_energy(), pose5b->energies().total_energy(), std::max( std::abs( std::max(pose5a->energies().total_energy(), pose5b->energies().total_energy())/1000.0 ), 1e-12 ) );
						if ( std::abs( pose5a->energies().total_energy() - pose5b->energies().total_energy() ) > std::max( std::abs( std::max(pose5a->energies().total_energy(), pose5b->energies().total_energy())/1000.0 ), 1e-12 ) ) {
							TR << "Failure on phi=" << iphi << " psi=" << ipsi << " omega=" << iphi << ".  Dumping..." << std::endl;
							char outstr[256];
							sprintf( outstr, "post_failure_A_%03i_%03i.pdb", iphi, ipsi );
							pose5a->dump_pdb(std::string(outstr));
							sprintf( outstr, "post_failure_B_%03i_%03i.pdb", iphi, ipsi );
							core::pose::Pose pose5bprime( *(pose5b->clone()) );
							for ( core::Size ir(1), irmax(pose5bprime.total_residue()); ir<=irmax; ++ir ) {
								for ( core::Size ia(1), iamax( pose5bprime.residue_type(ir).natoms()); ia<=iamax; ++ia ) {
									core::id::AtomID curat( ia, ir );
									numeric::xyzVector< core::Real > flipped_pos( pose5bprime.xyz(curat));
									flipped_pos.z( -1.0 * flipped_pos.z() );
									pose5bprime.set_xyz( curat, flipped_pos );
								}
							}
							pose5bprime.update_residue_neighbors();
							pose5bprime.dump_pdb(std::string(outstr));
						}
					}
					for ( core::Size ir=1, irmax=pose5a->size(); ir<=irmax; ++ir ) {
						TR << "iphi\t" << iphi << "\tipsi\t" << ipsi << std::endl;
						TR << "phi5a\t" << pose5a->phi(ir) << "\tphi5b\t" << pose5b->phi(ir) << std::endl;
						TR << "psi5a\t" << pose5a->psi(ir) << "\tpsi5b\t" << pose5b->psi(ir) << std::endl;
						TR << "omega5a\t" << pose5a->omega(ir) << "\tomega5b\t" << pose5b->omega(ir) << std::endl;
						if ( !skip_pose3_tests ) {
							TS_ASSERT( within_thresh( pose5a->phi(ir), -1.0*pose5b->phi(ir), 0.5 ) );
							TS_ASSERT( within_thresh( pose5a->psi(ir), -1.0*pose5b->psi(ir), 0.5 ) );
							TS_ASSERT( within_thresh( pose5a->omega(ir), -1.0*pose5b->omega(ir), 0.5 ) );
						}
					}


				} //if cartesian

			}
		}
	}

	/// @brief Tests symmetric minimization of A SINGLE glycine with the cart_bonded scorefunction.
	void test_single_gly_min_cart_bonded() {
		//Set up the scorefunction:
		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::cart_bonded, 1.0 );
		TR << "Testing cart_bonded score term on a single glycine." << std::endl;

		for ( int iphi=-180; iphi<=180; iphi+=60 ) {
			for ( int ipsi=-180; ipsi<=180; ipsi+=60 ) {
				core::pose::PoseOP pose( new core::pose::Pose() );
				core::pose::make_pose_from_sequence(*pose, "G", "fa_standard", false);
				core::pose::PoseOP pose2( pose->clone() );
				pose->set_phi(1, static_cast<core::Real>(iphi));
				pose->set_psi(1, static_cast<core::Real>(ipsi));
				pose2->set_phi(1, static_cast<core::Real>(-1*iphi));
				pose2->set_psi(1, static_cast<core::Real>(-1*ipsi));
				(*sfxn)(*pose);
				do_minimization(*pose, sfxn, true);
				(*sfxn)(*pose2);
				do_minimization(*pose2, sfxn, true);
				TS_ASSERT_DELTA( pose->energies().total_energy(), pose2->energies().total_energy(), std::max( std::abs( std::max(pose->energies().total_energy(), pose2->energies().total_energy())/1000.0 ), 1e-12 ) );
			}
		}
	}

	/// @brief Tests symmetric minimization of TWO glycines with the cart_bonded scorefunction.
	void test_two_gly_min_cart_bonded() {
		//Set up the scorefunction:
		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::cart_bonded, 1.0 );
		TR << "Testing cart_bonded score term on two glycines." << std::endl;

		for ( int iphi=-180; iphi<=180; iphi+=60 ) {
			for ( int ipsi=-180; ipsi<=180; ipsi+=60 ) {
				TR << "iphi=" << iphi << "\tipsi=" << ipsi << std::endl;
				core::pose::PoseOP pose( new core::pose::Pose() );
				core::pose::make_pose_from_sequence(*pose, "GG", "fa_standard", false);
				core::pose::PoseOP pose2( pose->clone() );
				pose->set_phi(1, static_cast<core::Real>(iphi));
				pose->set_psi(1, static_cast<core::Real>(ipsi));
				pose->set_omega(1, 175.0);
				pose->set_phi(2, static_cast<core::Real>(iphi));
				pose->set_psi(2, static_cast<core::Real>(ipsi));
				pose2->set_phi(1, -1.0*static_cast<core::Real>(iphi));
				pose2->set_psi(1, -1.0*static_cast<core::Real>(ipsi));
				pose2->set_omega(1, -175.0);
				pose2->set_phi(2, -1.0*static_cast<core::Real>(iphi));
				pose2->set_psi(2, -1.0*static_cast<core::Real>(ipsi));
				(*sfxn)(*pose);
				do_minimization(*pose, sfxn, true);
				(*sfxn)(*pose2);
				do_minimization(*pose2, sfxn, true);
				TS_ASSERT_DELTA( pose->energies().total_energy(), pose2->energies().total_energy(), std::max( std::abs( std::max(pose->energies().total_energy(), pose2->energies().total_energy())/1000.0 ), 1e-12 ) );
				TS_ASSERT_DELTA( pose->psi(1), -1.0*pose2->psi(1), 0.01 );
				TS_ASSERT_DELTA( pose->omega(1), -1.0*pose2->omega(1), 0.01 );
				TS_ASSERT_DELTA( pose->phi(2), -1.0*pose2->phi(2), 0.01 );
			}
		}
	}

	/// @brief Tests symmetric scoring of glycine with the cart_bonded scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_gly_min_cart_bonded() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::cart_bonded, 1.0 );
		TR << "Testing cart_bonded score term." << std::endl;
		repeat_structure_test(scorefxn, true);
		return;
	}

	/// @brief Tests symmetric scoring of glycine with the fa_atr scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_gly_min_fa_atr() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_atr, 1.0 );
		TR << "Testing fa_atr score term." << std::endl;
		repeat_structure_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring of glycine with the fa_rep scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_gly_min_fa_rep() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_rep, 1.0 );
		TR << "Testing fa_rep score term." << std::endl;
		repeat_structure_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring of glycine with the fa_intra_rep scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_gly_min_fa_intra_rep() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_intra_rep, 1.0 );
		TR << "Testing fa_intra_rep score term." << std::endl;
		repeat_structure_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring of glycine with the fa_sol scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_gly_min_fa_sol() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing fa_sol score term." << std::endl;
		repeat_structure_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring of glycine with the fa_elec scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_gly_min_fa_elec() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_elec, 1.0 );
		TR << "Testing fa_elec score term." << std::endl;
		repeat_structure_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring of glycine with the hbonds scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_gly_min_hbonds() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::hbond_sr_bb, 1.0 );
		scorefxn->set_weight( core::scoring::hbond_lr_bb, 1.0 );
		scorefxn->set_weight( core::scoring::hbond_sc, 1.0 );
		scorefxn->set_weight( core::scoring::hbond_bb_sc, 1.0 );
		TR << "Testing hbonds score terms." << std::endl;
		repeat_structure_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring of glycine with the fa_dun scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_gly_min_fa_dun() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_dun, 1.0 );
		TR << "Testing fa_dun score term." << std::endl;
		repeat_structure_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring of glycine with the omega scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_gly_min_omega() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::omega, 1.0 );
		TR << "Testing omega score term." << std::endl;
		repeat_structure_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring of glycine with the rama scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_gly_min_rama() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::rama, 1.0 );
		TR << "Testing rama score term." << std::endl;
		repeat_structure_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring of glycine with the rama_prepro scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_gly_min_rama_prepro() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::rama_prepro, 1.0 );
		TR << "Testing rama_prepro score term." << std::endl;
		repeat_structure_test(scorefxn, false);
		return;
	}


	/// @brief Tests symmetric scoring of glycine with the p_aa_pp scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_gly_min_p_aa_pp() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::p_aa_pp, 1.0 );
		TR << "Testing p_aa_pp score term." << std::endl;
		repeat_structure_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring of glycine with the full talaris2014 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_gly_min_talaris2014() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->add_weights_from_file("talaris2014.wts");
		TR << "Testing full talaris2014 score function." << std::endl;
		repeat_structure_test(scorefxn, false);
		return;
	}

	/// @brief Tests symmetric scoring of glycine with the full talaris2014_cart scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_gly_min_talaris2014_cart() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->add_weights_from_file("talaris2014_cart.wts");
		TR << "Testing full talaris2014_cart score function." << std::endl;
		repeat_structure_test(scorefxn, true);
		return;
	}

	/// @brief Tests symmetric scoring of glycine with the full default scorefunction, whatever that currently is.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_symm_gly_min_default_scorefxn() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
		TR << "Testing full default score function." << std::endl;
		repeat_structure_test(scorefxn, false);
		return;
	}

};
