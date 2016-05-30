// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/use_bicubic_interpolation.cxxtest.hh
/// @brief  Test the usage of the -corrections:score:use_bicubic_interpolation flag
/// @author Matthew O'Meara

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>


#include <core/types.hh>
#include <core/scoring/Ramachandran.hh>

// Unit headers
#include <core/scoring/EnergyMap.hh>

#include <test/UTracer.hh>

//Auto Headers
#include <platform/types.hh>

#include <core/scoring/ScoreType.hh>
#include <sstream>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.use_bicubic_interpolation.cxxtest");

class UseBicubicInterpolationTest : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_rama_maps() {
		do_test_rama_map("scoring/score_functions/rama/Rama_smooth_dyn.dat_ss_6.4");
		do_test_rama_map("scoring/score_functions/rama/Rama.10.2009.yfsong.dat");
	}

	void do_test_rama_map(
		std::string const & rama_map_filename,
		core::Real sample_step_size = 5
	) {
		TR
			<< "Testing bicubic interplotion on the rama map: '"
			<< rama_map_filename << "'" << std::endl;

		using core::Real;
		using core::Size;
		using namespace core::chemical;

		core::scoring::Ramachandran rama_term(rama_map_filename, true);
		Real phi_min(0), phi_max(360), psi_min(0), psi_max(360);
		for ( Size aa = 1; aa <= num_canonical_aas; aa++ ) {
			for ( Real phi=phi_min; phi <= phi_max; phi = phi + sample_step_size ) {
				for ( Real psi=psi_min; psi <= psi_max; psi = psi + sample_step_size ) {

					Real rama(0), drama_dphi(0), drama_dpsi(0);
					rama_term.eval_rama_score_residue(
						true, false,
						static_cast<AA>(aa), phi, psi, rama, drama_dphi, drama_dpsi);
					TS_ASSERT_EQUALS(rama, rama);
					TS_ASSERT_EQUALS(drama_dphi, drama_dphi);
					TS_ASSERT_EQUALS(drama_dpsi, drama_dpsi);
				}
				TR << std::endl;
			}
		}
	}

	void test_specific_parts_rama_map() {
		do_test_specific_parts_rama_map(
			"scoring/score_functions/rama/Rama.10.2009.yfsong.dat",
			core::chemical::aa_pro,
			289.51369082867473,
			158.66287604720782);
	}

	void do_test_specific_parts_rama_map(
		std::string const & rama_map_filename,
		core::chemical::AA aa,
		core::Real phi,
		core::Real psi
	) {
		TR
			<< "Testing bicubic interplotion on the rama map: '"
			<< rama_map_filename << "' " << std::endl;

		using core::Real;
		using core::Size;
		using namespace core::chemical;

		core::scoring::Ramachandran rama_term(rama_map_filename, true);

		Real rama, drama_dphi, drama_dpsi;
		rama_term.eval_rama_score_residue(
			true, false,
			aa, phi, psi, rama, drama_dphi, drama_dpsi);
		TS_ASSERT_EQUALS(rama, rama);
		TS_ASSERT_EQUALS(drama_dphi, drama_dphi);
		TS_ASSERT_EQUALS(drama_dpsi, drama_dpsi);
		TR
			<< "at (phi=" << phi << ", psi=" << psi << "): "
			<< rama << " " << drama_dphi << " " << drama_dpsi << std::endl;
	}

};
