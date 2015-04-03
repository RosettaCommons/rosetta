// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/conformation/ppo_torsion_bin.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/AA.hh>

using namespace std;

using core::Size;
using core::Real;
using core::pose::Pose;
using core::chemical::AA;
using core::scoring::Ramachandran;

AA amino_acids[] = { // {{{1
		core::chemical::aa_unk, // Start counting from 1.
		core::chemical::aa_ala,
		core::chemical::aa_cys,
		core::chemical::aa_asp,
		core::chemical::aa_glu,
		core::chemical::aa_phe,
		core::chemical::aa_gly,
		core::chemical::aa_his,
		core::chemical::aa_ile,
		core::chemical::aa_lys,
		core::chemical::aa_leu,
		core::chemical::aa_met,
		core::chemical::aa_asn,
		core::chemical::aa_pro,
		core::chemical::aa_gln,
		core::chemical::aa_arg,
		core::chemical::aa_ser,
		core::chemical::aa_thr,
		core::chemical::aa_val,
		core::chemical::aa_trp,
		core::chemical::aa_tyr };
// }}}1

// This is not a particularly strong set of tests, but hopefully it is enough
// to catch really bad mistakes.  Hopefully it will also serve as a good
// framework for building more tests in the future.

class RamachandranTest : public CxxTest::TestSuite {

public:

	// Create a mock pose and fill it with artificial phi/psi values.  This pose
	// is used by several of the tests to compare against.
	void setUp() {
		core_init();
		pose_ = core::pose::PoseOP( new core::pose::Pose );
		core::pose::make_pose_from_sequence(*pose_, "VVAV", "fa_standard", false);

		pose_->set_phi(1, 296); pose_->set_psi(1, 319);		// alpha helix
		pose_->set_phi(2, 235); pose_->set_psi(2, 138);		// beta strand
		pose_->set_phi(3,  55); pose_->set_psi(3,  42);		// left-handed helix
		pose_->set_phi(4, 125); pose_->set_psi(4, 100);		// forbidden
	}

	void tearDown() {
		pose_.reset();
	}

	// Test the score function on a handful of points, as defined in setup().
	void test_eval_rama_score_residue() {
		Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();
		Real expected[] = { -0.2578, -0.9390, 0.4680, 4.9683};

		for (Size i = 1; i <= pose_->total_residue(); i++) {
			Real observed = rama.eval_rama_score_residue(pose_->residue(i));
			TS_ASSERT_DELTA(observed, expected[i-1], 1e-4);
		}
	}

	// Just sample a few phi/psi pairs to make sure nothing is too broken.
	void test_random_phipsi_from_rama() {
		Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();
		Real phi, psi;

		for (Size i = 1; i <= 20; i++) {
			rama.random_phipsi_from_rama(amino_acids[i], phi, psi);
			//std::cout << "phipsi: " << phi << " " << psi << std::endl;
		}
	}

	// Just sample a few phi/psi pairs to make sure nothing is too broken.
	void test_uniform_phipsi_from_allowed_rama() {
		Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();
		Real phi, psi;

		for (Size i = 1; i <= 20; i++) {
			rama.uniform_phipsi_from_allowed_rama(amino_acids[i], phi, psi);
		}
	}

	// Make sure allowed and forbidden points are properly discriminated.
	void test_phipsi_in_allowed_rama() {
		Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();
		bool expected[] = {true, true, true, false};

		for (Size i = 1; i <= pose_->total_residue(); i++) {
			bool is_allowed = rama.phipsi_in_allowed_rama(
					pose_->aa(i), pose_->phi(i), pose_->psi(i));
			bool is_forbidden = rama.phipsi_in_forbidden_rama(
					pose_->aa(i), pose_->phi(i), pose_->psi(i));

			TS_ASSERT_EQUALS(is_allowed, expected[i-1]);
			TS_ASSERT_EQUALS(is_forbidden, not expected[i-1]);
		}
	}

	void test_rama_cdf_creation() {
		// first we're going to reconstruct the CDF by walking along the
		// phi/psi boundaries, and then we're going to make sure that the
		// cdf has been calculated correctly.
		Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();
		ObjexxFCL::FArray2D< Real > phe_prob( rama.n_phi_bins(), rama.n_psi_bins(), 0.0 );

		core::Real binsize = 360 / rama.n_phi_bins();

		core::Real probsum = 0;
		core::Real const minprob = rama.minimum_sampling_probability();
		for ( core::Size ii = 0; ii < rama.n_phi_bins(); ++ii ) {
			for ( core::Size jj = 0; jj < rama.n_psi_bins(); ++jj ) {
				phe_prob( ii+1, jj+1 ) = rama.rama_probability( core::chemical::aa_phe, ii*binsize, jj*binsize );
				if ( phe_prob( ii+1, jj+1 ) >= minprob ) {
					probsum += phe_prob( ii+1, jj+1 );
				}
			}
		}
		Real inv_probsum = 1 / probsum;

		utility::vector1< core::Real > const & phe_cdf = rama.cdf_for_aa( core::chemical::aa_phe );
		core::Size nphipsi_bins = rama.n_phi_bins() * rama.n_psi_bins();
		core::Size count = 1;
		for ( core::Size ii = 1; ii <= rama.n_phi_bins(); ++ii ) {
			for ( core::Size jj = 1; jj <= rama.n_psi_bins(); ++jj ) {
				if ( phe_prob( ii, jj ) >= minprob ) {
					if ( count < nphipsi_bins ) {
						TS_ASSERT( phe_cdf[ count ] != phe_cdf[ count+1 ] );
						//std::cout << " " << count << " " << phe_cdf[ count ] + phe_prob( ii, jj ) * inv_probsum << " " << phe_cdf[ count+1 ] << std::endl;
						TS_ASSERT_DELTA( phe_cdf[ count ] + phe_prob( ii, jj ) * inv_probsum, phe_cdf[ count+1 ], 1e-6 );
					} else {
						TS_ASSERT_DELTA( phe_cdf[ count ] + phe_prob( ii,jj ) * inv_probsum, 1.0, 1e-6 );
					}
				} else {
					if ( count < nphipsi_bins ) {
						//std::cout << " phe_prob( ii, jj ) < minprob: " << phe_cdf[ count ] << " " << phe_cdf[ count+1 ] << std::endl;
						TS_ASSERT( phe_cdf[ count ] == phe_cdf[ count+1 ] );
					} else {
						TS_ASSERT_DELTA( phe_cdf[ count ], 1.0, 1e-6 );
					}
				}
				++count;
			}
		}
	}

	void test_rama_cdf_creation_for_ppo_torbin_A() {
		// first we're going to reconstruct the CDF by walking along the
		// phi/psi boundaries, and then we're going to make sure that the
		// cdf has been calculated correctly.
		Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();
		ObjexxFCL::FArray2D< Real > phe_prob( rama.n_phi_bins(), rama.n_psi_bins(), 0.0 );

		core::Real binsize = 360 / rama.n_phi_bins();

		core::Real probsum = 0;
		core::Real const minprob = rama.minimum_sampling_probability();
		core::Size count_good = 0;
		for ( core::Size ii = 0; ii < rama.n_phi_bins(); ++ii ) {
			for ( core::Size jj = 0; jj < rama.n_psi_bins(); ++jj ) {
				phe_prob( ii+1, jj+1 ) = rama.rama_probability( core::chemical::aa_phe, ii*binsize, jj*binsize );
				core::conformation::ppo_torsion_bin iijjbin = core::conformation::get_torsion_bin( ii*binsize, jj*binsize, 180 );
				if ( phe_prob( ii+1, jj+1 ) >= minprob && iijjbin == core::conformation::ppo_torbin_A ) {
					probsum += phe_prob( ii+1, jj+1 );
					++count_good;
					// std::cout << "  " << ii*binsize << " " << jj*binsize << " " << phe_prob(ii+1, jj+1) << " " << core::conformation::map_torsion_bin_to_char( iijjbin ) << std::endl;
				}
			}
		}

		//std::cout << "count_good: " << count_good << std::endl;
		TS_ASSERT( count_good > 0 );

		Real inv_probsum = 1 / probsum;
		//std::cout << "inv_probsum: " << inv_probsum << std::endl;

		utility::vector1< core::Real > const & phe_cdf = rama.cdf_for_aa_for_torsion_bin( core::chemical::aa_phe, core::conformation::ppo_torbin_A );
		core::Size nphipsi_bins = rama.n_phi_bins() * rama.n_psi_bins();
		core::Size count = 1;
		for ( core::Size ii = 1; ii <= rama.n_phi_bins(); ++ii ) {
			for ( core::Size jj = 1; jj <= rama.n_psi_bins(); ++jj ) {
				core::conformation::ppo_torsion_bin iijjbin = core::conformation::get_torsion_bin( (ii-1)*binsize, (jj-1)*binsize, 180 );

				if ( phe_prob( ii, jj ) >= minprob && iijjbin == core::conformation::ppo_torbin_A ) {
					if ( count < nphipsi_bins ) {
						TS_ASSERT( phe_cdf[ count ] != phe_cdf[ count+1 ] );
						//std::cout << " " << count << " " << phe_cdf[ count ] + phe_prob( ii, jj ) * inv_probsum << " " << phe_cdf[ count+1 ] << std::endl;
						TS_ASSERT_DELTA( phe_cdf[ count ] + phe_prob( ii, jj ) * inv_probsum, phe_cdf[ count+1 ], 1e-6 );
					} else {
						TS_ASSERT_DELTA( phe_cdf[ count ] + phe_prob( ii,jj ) * inv_probsum, 1.0, 1e-6 );
					}
				} else {
					if ( count < nphipsi_bins ) {
						//std::cout << " phe_prob( ii, jj ) < minprob or iijjbin ( " << core::conformation::map_torsion_bin_to_char( iijjbin ) << ") != ppo_torbin_A: " << phe_cdf[ count ] << " " << phe_cdf[ count+1 ] << std::endl;
						TS_ASSERT( phe_cdf[ count ] == phe_cdf[ count+1 ] );
					} else {
						TS_ASSERT_DELTA( phe_cdf[ count ], 1.0, 1e-6 );
					}
				}
				++count;
			}
		}
		//std::cout << "BLAH!" << std::endl;
	}

	void test_random_phipsi_from_rama_by_torsion_bin() {
		Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();
		core::Real binsize = 360 / rama.n_phi_bins();
		for ( core::Size ii = 1; ii < 100; ++ii ) {
			Real phi, psi;
			rama.random_phipsi_from_rama_by_torsion_bin( core::chemical::aa_phe, phi, psi, core::conformation::ppo_torbin_A );
			Real round_down_phi = floor( phi / binsize ) * binsize;
			Real round_down_psi = floor( psi / binsize ) * binsize;
			core::conformation::ppo_torsion_bin torbin =
				core::conformation::get_torsion_bin( round_down_phi, round_down_psi );
			// std::cout << " phi: " << phi << " psi: " << psi << " round phi: " << round_down_phi << " round psi " << round_down_psi;
			// std::cout << " torbin: " << core::conformation::map_torsion_bin_to_char( torbin ) << std::endl;
			TS_ASSERT_EQUALS( torbin, core::conformation::ppo_torbin_A );
			if ( torbin != core::conformation::ppo_torbin_A ) {
				core::conformation::ppo_torsion_bin torbin2 = core::conformation::get_torsion_bin( round_down_phi, round_down_psi-1 );
				std::cout << " previous psi: " << core::conformation::map_torsion_bin_to_char( torbin2 ) << std::endl;
			}
		}
	}


private:
	core::pose::PoseOP pose_;

};
