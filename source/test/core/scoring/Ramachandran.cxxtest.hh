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
#include <core/scoring/Ramachandran.hh>

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
		core::pose::make_pose_from_sequence(pose, "VVVV", "fa_standard", false);

		pose.set_phi(1, 296); pose.set_psi(1, 319);		// alpha helix
		pose.set_phi(2, 235); pose.set_psi(2, 138);		// beta strand
		pose.set_phi(3,  55); pose.set_psi(3,  42);		// left-handed helix
		pose.set_phi(4, 125); pose.set_psi(4, 100);		// forbidden
	}

	void tearDown() {}

	// Test the score function on a handful of points, as defined in setup().
	void test_eval_rama_score_residue() {
		Ramachandran rama;
		Real expected[] = {0, -0.2578, -0.9390, 2.7557, 4.9683};

		for (Size i = 1; i <= pose.total_residue(); i++) {
			Real observed = rama.eval_rama_score_residue(pose.residue(i));
			TS_ASSERT_DELTA(observed, expected[i], 1e-4);
		}
	}

	// Just sample a few phi/psi pairs to make sure nothing is too broken.
	void test_random_phipsi_from_rama() {
		Ramachandran rama;
		Real phi, psi;

		for (Size i = 1; i <= 20; i++) {
			rama.random_phipsi_from_rama(amino_acids[i], phi, psi);
		}
	}

	// Just sample a few phi/psi pairs to make sure nothing is too broken.
	void test_uniform_phipsi_from_allowed_rama() {
		Ramachandran rama;
		Real phi, psi;

		for (Size i = 1; i <= 20; i++) {
			rama.uniform_phipsi_from_allowed_rama(amino_acids[i], phi, psi);
		}
	}

	// Make sure allowed and forbidden points are properly discriminated.
	void test_phipsi_in_allowed_rama() {
		Ramachandran rama;
		bool expected[] = {true, true, true, false};

		for (Size i = 1; i <= pose.total_residue(); i++) {
			bool is_allowed = rama.phipsi_in_allowed_rama(
					pose.aa(i), pose.phi(i), pose.psi(i));
			bool is_forbidden = rama.phipsi_in_forbidden_rama(
					pose.aa(i), pose.phi(i), pose.psi(i));

			TS_ASSERT_EQUALS(is_allowed, expected[i]);
			TS_ASSERT_EQUALS(is_forbidden, not expected[i]);
		}
	}

private:
	Pose pose;

};
