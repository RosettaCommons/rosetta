// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/conformation/ppo_torsion_bin.hh>
#include <core/scoring/Ramachandran2B.hh>
#include <core/scoring/ScoringManager.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/AA.hh>

using namespace std;

class Ramachandran2BTests : public CxxTest::TestSuite {

public:

	// Create a mock pose and fill it with artificial phi/psi values.  This pose
	// is used by several of the tests to compare against.
	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_rama2b_left_cdf_creation() {
		// first we're going to reconstruct the CDF by walking along the
		// phi/psi boundaries, and then we're going to make sure that the
		// cdf has been calculated correctly.
		core::scoring::Ramachandran2B const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran2B();
		ObjexxFCL::FArray2D< core::Real > asp_phe_prob( rama.n_phi_bins(), rama.n_psi_bins(), 0.0 );

		core::Real binsize = 360 / rama.n_phi_bins();

		core::Real probsum = 0;
		core::Real const minprob = rama.minimum_sampling_probability();
		for ( core::Size ii = 0; ii < rama.n_phi_bins(); ++ii ) {
			for ( core::Size jj = 0; jj < rama.n_psi_bins(); ++jj ) {
				asp_phe_prob( ii+1, jj+1 ) = rama.rama_bin_probability_left( core::chemical::aa_asp, core::chemical::aa_phe, ii*binsize, jj*binsize );
				if ( asp_phe_prob( ii+1, jj+1 ) >= minprob ) {
					probsum += asp_phe_prob( ii+1, jj+1 );
				}
			}
		}
		core::Real inv_probsum = 1 / probsum;

		utility::vector1< core::Real > const & asp_phe_cdf = rama.left_cdf( core::chemical::aa_asp, core::chemical::aa_phe );
		core::Size nphipsi_bins = rama.n_phi_bins() * rama.n_psi_bins();
		core::Size count = 1;
		for ( core::Size ii = 1; ii <= rama.n_phi_bins(); ++ii ) {
			for ( core::Size jj = 1; jj <= rama.n_psi_bins(); ++jj ) {
				if ( asp_phe_prob( ii, jj ) >= minprob ) {
					if ( count < nphipsi_bins ) {
						TS_ASSERT( asp_phe_cdf[ count ] != asp_phe_cdf[ count+1 ] );
						//std::cout << " " << count << " " << phe_cdf[ count ] + phe_prob( ii, jj ) * inv_probsum << " " << phe_cdf[ count+1 ] << std::endl;
						TS_ASSERT_DELTA( asp_phe_cdf[ count ] + asp_phe_prob( ii, jj ) * inv_probsum, asp_phe_cdf[ count+1 ], 1e-6 );
					} else {
						TS_ASSERT_DELTA( asp_phe_cdf[ count ] + asp_phe_prob( ii,jj ) * inv_probsum, 1.0, 1e-6 );
					}
				} else {
					if ( count < nphipsi_bins ) {
						//std::cout << " asp_phe_prob( ii, jj ) < minprob: " << asp_phe_cdf[ count ] << " " << asp_phe_cdf[ count+1 ] << std::endl;
						TS_ASSERT( asp_phe_cdf[ count ] == asp_phe_cdf[ count+1 ] );
					} else {
						TS_ASSERT_DELTA( asp_phe_cdf[ count ], 1.0, 1e-6 );
					}
				}
				++count;
			}
		}
	}

	void test_rama2b_right_cdf_creation() {
		// first we're going to reconstruct the CDF by walking along the
		// phi/psi boundaries, and then we're going to make sure that the
		// cdf has been calculated correctly.
		core::scoring::Ramachandran2B const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran2B();
		ObjexxFCL::FArray2D< core::Real > phe_asp_prob( rama.n_phi_bins(), rama.n_psi_bins(), 0.0 );

		core::Real binsize = 360 / rama.n_phi_bins();

		core::Real probsum = 0;
		core::Real const minprob = rama.minimum_sampling_probability();
		for ( core::Size ii = 0; ii < rama.n_phi_bins(); ++ii ) {
			for ( core::Size jj = 0; jj < rama.n_psi_bins(); ++jj ) {
				phe_asp_prob( ii+1, jj+1 ) = rama.rama_bin_probability_right( core::chemical::aa_phe, core::chemical::aa_asp, ii*binsize, jj*binsize );
				if ( phe_asp_prob( ii+1, jj+1 ) >= minprob ) {
					probsum += phe_asp_prob( ii+1, jj+1 );
				}
			}
		}
		core::Real inv_probsum = 1 / probsum;

		utility::vector1< core::Real > const & phe_asp_cdf = rama.right_cdf( core::chemical::aa_phe, core::chemical::aa_asp );
		core::Size nphipsi_bins = rama.n_phi_bins() * rama.n_psi_bins();
		core::Size count = 1;
		for ( core::Size ii = 1; ii <= rama.n_phi_bins(); ++ii ) {
			for ( core::Size jj = 1; jj <= rama.n_psi_bins(); ++jj ) {
				if ( phe_asp_prob( ii, jj ) >= minprob ) {
					if ( count < nphipsi_bins ) {
						TS_ASSERT( phe_asp_cdf[ count ] != phe_asp_cdf[ count+1 ] );
						//std::cout << " " << count << " " << phe_cdf[ count ] + phe_prob( ii, jj ) * inv_probsum << " " << phe_cdf[ count+1 ] << std::endl;
						TS_ASSERT_DELTA( phe_asp_cdf[ count ] + phe_asp_prob( ii, jj ) * inv_probsum, phe_asp_cdf[ count+1 ], 1e-6 );
					} else {
						TS_ASSERT_DELTA( phe_asp_cdf[ count ] + phe_asp_prob( ii,jj ) * inv_probsum, 1.0, 1e-6 );
					}
				} else {
					if ( count < nphipsi_bins ) {
						//std::cout << " phe_asp_prob( ii, jj ) < minprob: " << phe_asp_cdf[ count ] << " " << phe_asp_cdf[ count+1 ] << std::endl;
						TS_ASSERT( phe_asp_cdf[ count ] == phe_asp_cdf[ count+1 ] );
					} else {
						TS_ASSERT_DELTA( phe_asp_cdf[ count ], 1.0, 1e-6 );
					}
				}
				++count;
			}
		}
	}

	void test_rama2b_left_cdf_creation_for_ppo_torbin_A() {
		// first we're going to reconstruct the CDF by walking along the
		// phi/psi boundaries, and then we're going to make sure that the
		// cdf has been calculated correctly.
		core::scoring::Ramachandran2B const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran2B();
		ObjexxFCL::FArray2D< core::Real > asp_phe_prob( rama.n_phi_bins(), rama.n_psi_bins(), 0.0 );

		core::Real binsize = 360 / rama.n_phi_bins();

		core::Real probsum = 0;
		core::Real const minprob = rama.minimum_sampling_probability();
		core::Size count_good = 0;
		for ( core::Size ii = 0; ii < rama.n_phi_bins(); ++ii ) {
			for ( core::Size jj = 0; jj < rama.n_psi_bins(); ++jj ) {
				asp_phe_prob( ii+1, jj+1 ) = rama.rama_bin_probability_left( core::chemical::aa_asp, core::chemical::aa_phe, ii*binsize, jj*binsize );
				core::conformation::ppo_torsion_bin iijjbin = core::conformation::get_torsion_bin( ii*binsize, jj*binsize, 180 );
				if ( asp_phe_prob( ii+1, jj+1 ) >= minprob && iijjbin == core::conformation::ppo_torbin_A ) {
					probsum += asp_phe_prob( ii+1, jj+1 );
					++count_good;
					// std::cout << "  " << ii*binsize << " " << jj*binsize << " " << asp_phe_prob(ii+1, jj+1) << " " << core::conformation::map_torsion_bin_to_char( iijjbin ) << std::endl;
				}
			}
		}

		//std::cout << "count_good: " << count_good << std::endl;
		TS_ASSERT( count_good > 0 );

		core::Real inv_probsum = 1 / probsum;
		//std::cout << "inv_probsum: " << inv_probsum << std::endl;

		utility::vector1< core::Real > const & asp_phe_cdf =
			rama.left_cdf_for_torsion_bin( core::chemical::aa_asp, core::chemical::aa_phe, core::conformation::ppo_torbin_A );
		core::Size nphipsi_bins = rama.n_phi_bins() * rama.n_psi_bins();
		core::Size count = 1;
		for ( core::Size ii = 1; ii <= rama.n_phi_bins(); ++ii ) {
			for ( core::Size jj = 1; jj <= rama.n_psi_bins(); ++jj ) {
				core::conformation::ppo_torsion_bin iijjbin = core::conformation::get_torsion_bin( (ii-1)*binsize, (jj-1)*binsize, 180 );

				if ( asp_phe_prob( ii, jj ) >= minprob && iijjbin == core::conformation::ppo_torbin_A ) {
					if ( count < nphipsi_bins ) {
						TS_ASSERT( asp_phe_cdf[ count ] != asp_phe_cdf[ count+1 ] );
						//std::cout << " " << count << " " << asp_phe_cdf[ count ] + asp_phe_prob( ii, jj ) * inv_probsum << " " << asp_phe_cdf[ count+1 ] << std::endl;
						TS_ASSERT_DELTA( asp_phe_cdf[ count ] + asp_phe_prob( ii, jj ) * inv_probsum, asp_phe_cdf[ count+1 ], 1e-6 );
					} else {
						TS_ASSERT_DELTA( asp_phe_cdf[ count ] + asp_phe_prob( ii,jj ) * inv_probsum, 1.0, 1e-6 );
					}
				} else {
					if ( count < nphipsi_bins ) {
						//std::cout << " asp_phe_prob( ii, jj ) < minprob or iijjbin ( " << core::conformation::map_torsion_bin_to_char( iijjbin ) << ") != ppo_torbin_A: " << asp_phe_cdf[ count ] << " " << asp_phe_cdf[ count+1 ] << std::endl;
						TS_ASSERT( asp_phe_cdf[ count ] == asp_phe_cdf[ count+1 ] );
					} else {
						TS_ASSERT_DELTA( asp_phe_cdf[ count ], 1.0, 1e-6 );
					}
				}
				++count;
			}
		}
		//std::cout << "BLAH!" << std::endl;
	}

	void test_rama2b_right_cdf_creation_for_ppo_torbin_A() {
		// first we're going to reconstruct the CDF by walking along the
		// phi/psi boundaries, and then we're going to make sure that the
		// cdf has been calculated correctly.
		core::scoring::Ramachandran2B const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran2B();
		ObjexxFCL::FArray2D< core::Real > asp_phe_prob( rama.n_phi_bins(), rama.n_psi_bins(), 0.0 );

		core::Real binsize = 360 / rama.n_phi_bins();

		core::Real probsum = 0;
		core::Real const minprob = rama.minimum_sampling_probability();
		core::Size count_good = 0;
		for ( core::Size ii = 0; ii < rama.n_phi_bins(); ++ii ) {
			for ( core::Size jj = 0; jj < rama.n_psi_bins(); ++jj ) {
				asp_phe_prob( ii+1, jj+1 ) = rama.rama_bin_probability_right( core::chemical::aa_asp, core::chemical::aa_phe, ii*binsize, jj*binsize );
				core::conformation::ppo_torsion_bin iijjbin = core::conformation::get_torsion_bin( ii*binsize, jj*binsize, 180 );
				if ( asp_phe_prob( ii+1, jj+1 ) >= minprob && iijjbin == core::conformation::ppo_torbin_A ) {
					probsum += asp_phe_prob( ii+1, jj+1 );
					++count_good;
					// std::cout << "  " << ii*binsize << " " << jj*binsize << " " << asp_phe_prob(ii+1, jj+1) << " " << core::conformation::map_torsion_bin_to_char( iijjbin ) << std::endl;
				}
			}
		}

		//std::cout << "count_good: " << count_good << std::endl;
		TS_ASSERT( count_good > 0 );

		core::Real inv_probsum = 1 / probsum;
		//std::cout << "inv_probsum: " << inv_probsum << std::endl;

		utility::vector1< core::Real > const & asp_phe_cdf =
			rama.right_cdf_for_torsion_bin( core::chemical::aa_asp, core::chemical::aa_phe, core::conformation::ppo_torbin_A );
		core::Size nphipsi_bins = rama.n_phi_bins() * rama.n_psi_bins();
		core::Size count = 1;
		for ( core::Size ii = 1; ii <= rama.n_phi_bins(); ++ii ) {
			for ( core::Size jj = 1; jj <= rama.n_psi_bins(); ++jj ) {
				core::conformation::ppo_torsion_bin iijjbin = core::conformation::get_torsion_bin( (ii-1)*binsize, (jj-1)*binsize, 180 );

				if ( asp_phe_prob( ii, jj ) >= minprob && iijjbin == core::conformation::ppo_torbin_A ) {
					if ( count < nphipsi_bins ) {
						TS_ASSERT( asp_phe_cdf[ count ] != asp_phe_cdf[ count+1 ] );
						//std::cout << " " << count << " " << asp_phe_cdf[ count ] + asp_phe_prob( ii, jj ) * inv_probsum << " " << asp_phe_cdf[ count+1 ] << std::endl;
						TS_ASSERT_DELTA( asp_phe_cdf[ count ] + asp_phe_prob( ii, jj ) * inv_probsum, asp_phe_cdf[ count+1 ], 1e-6 );
					} else {
						TS_ASSERT_DELTA( asp_phe_cdf[ count ] + asp_phe_prob( ii,jj ) * inv_probsum, 1.0, 1e-6 );
					}
				} else {
					if ( count < nphipsi_bins ) {
						//std::cout << " asp_phe_prob( ii, jj ) < minprob or iijjbin ( " << core::conformation::map_torsion_bin_to_char( iijjbin ) << ") != ppo_torbin_A: " << asp_phe_cdf[ count ] << " " << asp_phe_cdf[ count+1 ] << std::endl;
						TS_ASSERT( asp_phe_cdf[ count ] == asp_phe_cdf[ count+1 ] );
					} else {
						TS_ASSERT_DELTA( asp_phe_cdf[ count ], 1.0, 1e-6 );
					}
				}
				++count;
			}
		}
		//std::cout << "BLAH!" << std::endl;
	}

	void test_random_phipsi_from_rama_by_torsion_bin_left() {
		core::scoring::Ramachandran2B const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran2B();
		core::Real binsize = 360 / rama.n_phi_bins();
		for ( core::Size ii = 1; ii < 100; ++ii ) {
			core::Real phi, psi;
			rama.random_phipsi_from_rama_by_torsion_bin_left( core::chemical::aa_asp, core::chemical::aa_phe, phi, psi, core::conformation::ppo_torbin_A );
			core::Real round_down_phi = floor( phi / binsize ) * binsize;
			core::Real round_down_psi = floor( psi / binsize ) * binsize;
			core::conformation::ppo_torsion_bin torbin =
				core::conformation::get_torsion_bin( round_down_phi, round_down_psi );
			//std::cout << " phi: " << phi << " psi: " << psi << " round phi: " << round_down_phi << " round psi " << round_down_psi;
			//std::cout << " torbin: " << core::conformation::map_torsion_bin_to_char( torbin ) << std::endl;
			TS_ASSERT_EQUALS( torbin, core::conformation::ppo_torbin_A );
			if ( torbin != core::conformation::ppo_torbin_A ) {
				core::conformation::ppo_torsion_bin torbin2 = core::conformation::get_torsion_bin( round_down_phi, round_down_psi-1 );
				std::cout << " previous psi: " << core::conformation::map_torsion_bin_to_char( torbin2 ) << std::endl;
			}
		}
	}

	void test_random_phipsi_from_rama_by_torsion_bin_right() {
		core::scoring::Ramachandran2B const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran2B();
		core::Real binsize = 360 / rama.n_phi_bins();
		for ( core::Size ii = 1; ii < 100; ++ii ) {
			core::Real phi, psi;
			rama.random_phipsi_from_rama_by_torsion_bin_right( core::chemical::aa_asp, core::chemical::aa_phe, phi, psi, core::conformation::ppo_torbin_A );
			core::Real round_down_phi = floor( phi / binsize ) * binsize;
			core::Real round_down_psi = floor( psi / binsize ) * binsize;
			core::conformation::ppo_torsion_bin torbin =
				core::conformation::get_torsion_bin( round_down_phi, round_down_psi );
			//std::cout << " phi: " << phi << " psi: " << psi << " round phi: " << round_down_phi << " round psi " << round_down_psi;
			//std::cout << " torbin: " << core::conformation::map_torsion_bin_to_char( torbin ) << std::endl;
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
