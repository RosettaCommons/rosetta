// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  numeric/crick_equations/BundleParams_derivatives.cxxtest.hh
/// @brief  Unit tests for analytical Crick equation derivatives.
/// @author Andy Watkins

#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <numeric/crick_equations/BundleParams.hh>
#include <numeric/crick_equations/BundleParams_derivatives.hh>

#include <cmath>

static basic::Tracer TR("CrickDerivativesTests");

class CrickDerivativesTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	numeric::xyzVector<numeric::Real> numerical_XYZ(
		numeric::Real t, numeric::Real r0, numeric::Real omega0, numeric::Real delta_omega0,
		numeric::Real r1, numeric::Real omega1, numeric::Real z1,
		numeric::Real delta_omega1, numeric::Real delta_z1
	) {
		bool failed = false;
		return numeric::crick_equations::XYZ_BUNDLE(
			t, r0, omega0, delta_omega0, r1, omega1, z1,
			delta_omega1, delta_z1, 1.0 /*epsilon*/, failed );
	}

	void check_derivative(
		std::string const & param_name,
		numeric::Real t, numeric::Real r0, numeric::Real omega0, numeric::Real delta_omega0,
		numeric::Real r1, numeric::Real omega1, numeric::Real z1,
		numeric::Real delta_omega1, numeric::Real delta_z1
	) {
		using numeric::Real;

		bool failed = false;
		numeric::crick_equations::CrickDerivatives derivs =
			numeric::crick_equations::compute_crick_derivatives(
				t, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1, failed );
		TS_ASSERT( !failed );
		if ( failed ) return;

		Real const eps = 1e-6;
		numeric::xyzVector<Real> analytical_deriv;
		numeric::xyzVector<Real> fwd, bwd;

		if ( param_name == "r0" ) {
			analytical_deriv = derivs.dXYZ_dr0;
			fwd = numerical_XYZ(t, r0+eps, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1);
			bwd = numerical_XYZ(t, r0-eps, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1);
		} else if ( param_name == "omega0" ) {
			analytical_deriv = derivs.dXYZ_domega0;
			fwd = numerical_XYZ(t, r0, omega0+eps, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1);
			bwd = numerical_XYZ(t, r0, omega0-eps, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1);
		} else if ( param_name == "delta_omega0" ) {
			analytical_deriv = derivs.dXYZ_ddelta_omega0;
			fwd = numerical_XYZ(t, r0, omega0, delta_omega0+eps, r1, omega1, z1, delta_omega1, delta_z1);
			bwd = numerical_XYZ(t, r0, omega0, delta_omega0-eps, r1, omega1, z1, delta_omega1, delta_z1);
		} else if ( param_name == "delta_omega1" ) {
			analytical_deriv = derivs.dXYZ_ddelta_omega1;
			fwd = numerical_XYZ(t, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1+eps, delta_z1);
			bwd = numerical_XYZ(t, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1-eps, delta_z1);
		} else if ( param_name == "delta_t" ) {
			analytical_deriv = derivs.dXYZ_ddelta_t;
			fwd = numerical_XYZ(t+eps, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1);
			bwd = numerical_XYZ(t-eps, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1);
		} else {
			TS_FAIL("Unknown parameter name: " + param_name);
			return;
		}

		numeric::xyzVector<Real> numerical_deriv = (fwd - bwd) / (2.0 * eps);

		Real const tol = 1e-4;
		TS_ASSERT_DELTA( analytical_deriv.x(), numerical_deriv.x(), tol );
		TS_ASSERT_DELTA( analytical_deriv.y(), numerical_deriv.y(), tol );
		TS_ASSERT_DELTA( analytical_deriv.z(), numerical_deriv.z(), tol );

		if ( std::abs(analytical_deriv.x() - numerical_deriv.x()) > tol ||
				std::abs(analytical_deriv.y() - numerical_deriv.y()) > tol ||
				std::abs(analytical_deriv.z() - numerical_deriv.z()) > tol ) {
			TR << "MISMATCH for " << param_name << " at t=" << t
				<< " r0=" << r0 << " omega0=" << omega0 << std::endl;
			TR << "  analytical: " << analytical_deriv.x() << " " << analytical_deriv.y() << " " << analytical_deriv.z() << std::endl;
			TR << "  numerical:  " << numerical_deriv.x() << " " << numerical_deriv.y() << " " << numerical_deriv.z() << std::endl;
		}
	}

	void test_derivatives_alpha_helix_params() {
		TR << "Testing Crick derivatives against numerical for alpha helix parameters." << std::endl;

		// Alpha helix-like parameters
		numeric::Real const r0 = 5.0;
		numeric::Real const omega0 = -0.045;  // radians per residue
		numeric::Real const delta_omega0 = 0.5;
		numeric::Real const r1 = 2.26;
		numeric::Real const omega1 = 1.72;
		numeric::Real const z1 = 1.51;
		numeric::Real const delta_omega1 = 0.0;
		numeric::Real const delta_z1 = 0.0;

		for ( numeric::Real t = -5.0; t <= 5.0; t += 2.5 ) {
			check_derivative("r0", t, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1);
			check_derivative("omega0", t, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1);
			check_derivative("delta_omega0", t, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1);
			check_derivative("delta_omega1", t, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1);
			check_derivative("delta_t", t, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1);
		}
	}

	void test_derivatives_with_nonzero_delta_z1() {
		TR << "Testing Crick derivatives with nonzero delta_z1." << std::endl;

		numeric::Real const r0 = 6.5;
		numeric::Real const omega0 = -0.045;
		numeric::Real const delta_omega0 = 1.2;
		numeric::Real const r1 = 1.5;
		numeric::Real const omega1 = 1.72;
		numeric::Real const z1 = 1.51;
		numeric::Real const delta_omega1 = 0.48;
		numeric::Real const delta_z1 = 1.05;

		for ( numeric::Real t = -3.0; t <= 3.0; t += 1.5 ) {
			check_derivative("r0", t, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1);
			check_derivative("omega0", t, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1);
			check_derivative("delta_omega0", t, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1);
			check_derivative("delta_omega1", t, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1);
			check_derivative("delta_t", t, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1);
		}
	}

};
