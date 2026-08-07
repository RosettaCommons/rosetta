// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/helical_bundle/PerturbBundleCorrelatedTests.cxxtest.hh
/// @brief  Unit tests for correlated perturbation in the PerturbBundle mover.
/// @author Andy Watkins

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>

// Project Headers
#include <protocols/helical_bundle/MakeBundle.hh>
#include <protocols/helical_bundle/PerturbBundle.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.hh>
#include <protocols/helical_bundle/parameters/BundleParametersSet.hh>
#include <protocols/helical_bundle/BundleParametrizationCalculator.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>

// Numeric Headers
#include <numeric/linear_algebra/cholesky_decomposition.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

#include <cmath>

static basic::Tracer TR("PerturbBundleCorrelatedTests");

class PerturbBundleCorrelatedTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	/// @brief Test that the Cholesky wrapper returns identity for identity input.
	void test_cholesky_identity() {
		TR << "Testing Cholesky decomposition of identity matrix." << std::endl;
		utility::vector1< utility::vector1< double > > identity( 3, utility::vector1< double >( 3, 0.0 ) );
		identity[1][1] = 1.0;
		identity[2][2] = 1.0;
		identity[3][3] = 1.0;

		utility::vector1< utility::vector1< double > > L = numeric::linear_algebra::cholesky_factor( identity );

		for ( platform::Size i = 1; i <= 3; ++i ) {
			for ( platform::Size j = 1; j <= 3; ++j ) {
				if ( i == j ) {
					TS_ASSERT_DELTA( L[i][j], 1.0, 1e-10 );
				} else {
					TS_ASSERT_DELTA( L[i][j], 0.0, 1e-10 );
				}
			}
		}
	}

	/// @brief Test Cholesky on a known 2x2 covariance matrix.
	void test_cholesky_2x2() {
		TR << "Testing Cholesky decomposition of 2x2 matrix." << std::endl;
		// Sigma = [[4, 2], [2, 3]]
		// L = [[2, 0], [1, sqrt(2)]]
		utility::vector1< utility::vector1< double > > cov( 2, utility::vector1< double >( 2, 0.0 ) );
		cov[1][1] = 4.0;
		cov[1][2] = 2.0;
		cov[2][1] = 2.0;
		cov[2][2] = 3.0;

		utility::vector1< utility::vector1< double > > L = numeric::linear_algebra::cholesky_factor( cov );

		TS_ASSERT_DELTA( L[1][1], 2.0, 1e-10 );
		TS_ASSERT_DELTA( L[1][2], 0.0, 1e-10 );
		TS_ASSERT_DELTA( L[2][1], 1.0, 1e-10 );
		TS_ASSERT_DELTA( L[2][2], std::sqrt(2.0), 1e-10 );
	}

	/// @brief Test that a non-positive-definite matrix throws.
	void test_cholesky_non_pd_throws() {
		TR << "Testing that non-PD matrix throws." << std::endl;
		utility::vector1< utility::vector1< double > > bad( 2, utility::vector1< double >( 2, 0.0 ) );
		bad[1][1] = 1.0;
		bad[1][2] = 2.0;
		bad[2][1] = 2.0;
		bad[2][2] = 1.0; // Not positive definite: det = 1 - 4 = -3

		TS_ASSERT_THROWS( numeric::linear_algebra::cholesky_factor( bad ), utility::excn::Exception & );
	}

	/// @brief Test perfect positive correlation (rho=1.0): two parameters should get identical deltas.
	void test_perfect_positive_correlation() {
		TR << "Testing perfect positive correlation (rho=1.0)." << std::endl;

		// Make a 2-helix bundle
		utility::tag::TagCOP maketag( tagptr_from_string(
			"<MakeBundle name=\"make\">\n<Helix r0=\"10\" />\n<Helix r0=\"10\" />\n</MakeBundle>") );
		// Perturb r0 of both helices with sigma=1.0 and correlation=1.0
		utility::tag::TagCOP perttag( tagptr_from_string(
			"<PerturbBundle name=\"pert\" default_perturbation_type=\"gaussian\">\n"
			"<Helix helix_index=\"1\" r0_perturbation=\"1.0\" />\n"
			"<Helix helix_index=\"2\" r0_perturbation=\"1.0\" />\n"
			"<Correlation helix1=\"1\" param1=\"r0\" helix2=\"2\" param2=\"r0\" correlation=\"1.0\" />\n"
			"</PerturbBundle>") );

		basic::datacache::DataMap dummy_data;

		for ( core::Size trial = 1; trial <= 10; ++trial ) {
			core::pose::Pose testpose;

			protocols::helical_bundle::MakeBundle makebundle;
			makebundle.parse_my_tag( maketag, dummy_data );
			makebundle.apply( testpose );

			core::Real r0_1_before, r0_2_before;
			{
				auto paramset = utility::pointer::dynamic_pointer_cast<
					protocols::helical_bundle::parameters::BundleParametersSet const >(
					testpose.conformation().parameters_set(1) );
				auto params1 = utility::pointer::dynamic_pointer_cast<
					protocols::helical_bundle::parameters::BundleParameters const >( paramset->parameters(1) );
				auto params2 = utility::pointer::dynamic_pointer_cast<
					protocols::helical_bundle::parameters::BundleParameters const >( paramset->parameters(2) );
				auto r0_p1 = utility::pointer::dynamic_pointer_cast<
					core::conformation::parametric::RealValuedParameter const >(
					params1->parameter_cop( protocols::helical_bundle::BPC_r0 ) );
				auto r0_p2 = utility::pointer::dynamic_pointer_cast<
					core::conformation::parametric::RealValuedParameter const >(
					params2->parameter_cop( protocols::helical_bundle::BPC_r0 ) );
				r0_1_before = r0_p1->value();
				r0_2_before = r0_p2->value();
			}

			protocols::helical_bundle::PerturbBundle pertbundle;
			pertbundle.parse_my_tag( perttag, dummy_data );
			pertbundle.apply( testpose );

			{
				auto paramset = utility::pointer::dynamic_pointer_cast<
					protocols::helical_bundle::parameters::BundleParametersSet const >(
					testpose.conformation().parameters_set(1) );
				auto params1 = utility::pointer::dynamic_pointer_cast<
					protocols::helical_bundle::parameters::BundleParameters const >( paramset->parameters(1) );
				auto params2 = utility::pointer::dynamic_pointer_cast<
					protocols::helical_bundle::parameters::BundleParameters const >( paramset->parameters(2) );
				auto r0_p1 = utility::pointer::dynamic_pointer_cast<
					core::conformation::parametric::RealValuedParameter const >(
					params1->parameter_cop( protocols::helical_bundle::BPC_r0 ) );
				auto r0_p2 = utility::pointer::dynamic_pointer_cast<
					core::conformation::parametric::RealValuedParameter const >(
					params2->parameter_cop( protocols::helical_bundle::BPC_r0 ) );

				core::Real delta1 = r0_p1->value() - r0_1_before;
				core::Real delta2 = r0_p2->value() - r0_2_before;

				// With rho=1.0 and equal sigmas, deltas should be identical
				TS_ASSERT_DELTA( delta1, delta2, 1e-10 );
			}
		}
	}

	/// @brief Test perfect negative correlation (rho=-1.0): deltas should have opposite signs.
	void test_perfect_negative_correlation() {
		TR << "Testing perfect negative correlation (rho=-1.0)." << std::endl;

		utility::tag::TagCOP maketag( tagptr_from_string(
			"<MakeBundle name=\"make\">\n<Helix r0=\"10\" />\n<Helix r0=\"10\" />\n</MakeBundle>") );
		utility::tag::TagCOP perttag( tagptr_from_string(
			"<PerturbBundle name=\"pert\" default_perturbation_type=\"gaussian\">\n"
			"<Helix helix_index=\"1\" r0_perturbation=\"1.0\" />\n"
			"<Helix helix_index=\"2\" r0_perturbation=\"1.0\" />\n"
			"<Correlation helix1=\"1\" param1=\"r0\" helix2=\"2\" param2=\"r0\" correlation=\"-1.0\" />\n"
			"</PerturbBundle>") );

		basic::datacache::DataMap dummy_data;

		for ( core::Size trial = 1; trial <= 10; ++trial ) {
			core::pose::Pose testpose;

			protocols::helical_bundle::MakeBundle makebundle;
			makebundle.parse_my_tag( maketag, dummy_data );
			makebundle.apply( testpose );

			core::Real r0_1_before, r0_2_before;
			{
				auto paramset = utility::pointer::dynamic_pointer_cast<
					protocols::helical_bundle::parameters::BundleParametersSet const >(
					testpose.conformation().parameters_set(1) );
				auto params1 = utility::pointer::dynamic_pointer_cast<
					protocols::helical_bundle::parameters::BundleParameters const >( paramset->parameters(1) );
				auto params2 = utility::pointer::dynamic_pointer_cast<
					protocols::helical_bundle::parameters::BundleParameters const >( paramset->parameters(2) );
				auto r0_p1 = utility::pointer::dynamic_pointer_cast<
					core::conformation::parametric::RealValuedParameter const >(
					params1->parameter_cop( protocols::helical_bundle::BPC_r0 ) );
				auto r0_p2 = utility::pointer::dynamic_pointer_cast<
					core::conformation::parametric::RealValuedParameter const >(
					params2->parameter_cop( protocols::helical_bundle::BPC_r0 ) );
				r0_1_before = r0_p1->value();
				r0_2_before = r0_p2->value();
			}

			protocols::helical_bundle::PerturbBundle pertbundle;
			pertbundle.parse_my_tag( perttag, dummy_data );
			pertbundle.apply( testpose );

			{
				auto paramset = utility::pointer::dynamic_pointer_cast<
					protocols::helical_bundle::parameters::BundleParametersSet const >(
					testpose.conformation().parameters_set(1) );
				auto params1 = utility::pointer::dynamic_pointer_cast<
					protocols::helical_bundle::parameters::BundleParameters const >( paramset->parameters(1) );
				auto params2 = utility::pointer::dynamic_pointer_cast<
					protocols::helical_bundle::parameters::BundleParameters const >( paramset->parameters(2) );
				auto r0_p1 = utility::pointer::dynamic_pointer_cast<
					core::conformation::parametric::RealValuedParameter const >(
					params1->parameter_cop( protocols::helical_bundle::BPC_r0 ) );
				auto r0_p2 = utility::pointer::dynamic_pointer_cast<
					core::conformation::parametric::RealValuedParameter const >(
					params2->parameter_cop( protocols::helical_bundle::BPC_r0 ) );

				core::Real delta1 = r0_p1->value() - r0_1_before;
				core::Real delta2 = r0_p2->value() - r0_2_before;

				// With rho=-1.0 and equal sigmas, deltas should be equal magnitude but opposite sign
				TS_ASSERT_DELTA( delta1, -delta2, 1e-10 );
			}
		}
	}

	/// @brief Test that invalid correlation (out of [-1,1]) throws during parsing.
	void test_invalid_correlation_throws() {
		TR << "Testing that invalid correlation coefficient throws." << std::endl;

		utility::tag::TagCOP perttag( tagptr_from_string(
			"<PerturbBundle name=\"pert\" default_perturbation_type=\"gaussian\">\n"
			"<Helix helix_index=\"1\" r0_perturbation=\"1.0\" />\n"
			"<Helix helix_index=\"2\" r0_perturbation=\"1.0\" />\n"
			"<Correlation helix1=\"1\" param1=\"r0\" helix2=\"2\" param2=\"r0\" correlation=\"1.5\" />\n"
			"</PerturbBundle>") );

		basic::datacache::DataMap dummy_data;
		protocols::helical_bundle::PerturbBundle pertbundle;
		TS_ASSERT_THROWS_ANYTHING( pertbundle.parse_my_tag( perttag, dummy_data ) );
	}

	/// @brief Test mixed mode: some parameters correlated, others perturbed independently.
	void test_mixed_correlated_and_independent() {
		TR << "Testing mixed correlated and independent perturbation." << std::endl;

		utility::tag::TagCOP maketag( tagptr_from_string(
			"<MakeBundle name=\"make\">\n<Helix r0=\"10\" />\n<Helix r0=\"10\" />\n</MakeBundle>") );
		// r0 on both helices correlated; omega0 on helix 1 independent
		utility::tag::TagCOP perttag( tagptr_from_string(
			"<PerturbBundle name=\"pert\" default_perturbation_type=\"gaussian\">\n"
			"<Helix helix_index=\"1\" r0_perturbation=\"1.0\" omega0_perturbation=\"0.01\" />\n"
			"<Helix helix_index=\"2\" r0_perturbation=\"1.0\" />\n"
			"<Correlation helix1=\"1\" param1=\"r0\" helix2=\"2\" param2=\"r0\" correlation=\"0.9\" />\n"
			"</PerturbBundle>") );

		basic::datacache::DataMap dummy_data;
		core::pose::Pose testpose;

		protocols::helical_bundle::MakeBundle makebundle;
		makebundle.parse_my_tag( maketag, dummy_data );
		makebundle.apply( testpose );

		core::Real omega0_before;
		{
			auto paramset = utility::pointer::dynamic_pointer_cast<
				protocols::helical_bundle::parameters::BundleParametersSet const >(
				testpose.conformation().parameters_set(1) );
			auto params1 = utility::pointer::dynamic_pointer_cast<
				protocols::helical_bundle::parameters::BundleParameters const >( paramset->parameters(1) );
			auto omega0_p1 = utility::pointer::dynamic_pointer_cast<
				core::conformation::parametric::RealValuedParameter const >(
				params1->parameter_cop( protocols::helical_bundle::BPC_omega0 ) );
			omega0_before = omega0_p1->value();
		}

		protocols::helical_bundle::PerturbBundle pertbundle;
		pertbundle.parse_my_tag( perttag, dummy_data );
		pertbundle.apply( testpose );

		{
			auto paramset = utility::pointer::dynamic_pointer_cast<
				protocols::helical_bundle::parameters::BundleParametersSet const >(
				testpose.conformation().parameters_set(1) );
			auto params1 = utility::pointer::dynamic_pointer_cast<
				protocols::helical_bundle::parameters::BundleParameters const >( paramset->parameters(1) );
			auto omega0_p1 = utility::pointer::dynamic_pointer_cast<
				core::conformation::parametric::RealValuedParameter const >(
				params1->parameter_cop( protocols::helical_bundle::BPC_omega0 ) );

			// omega0 should have been perturbed independently (not frozen)
			// This might occasionally fail if the perturbation is exactly 0, but with sigma=0.01 that's astronomically unlikely
			TS_ASSERT_DIFFERS( omega0_p1->value(), omega0_before );
		}
	}

	/// @brief Test that inconsistent correlations (non-PD matrix) throws during parsing.
	void test_inconsistent_correlations_throws() {
		TR << "Testing that mutually inconsistent correlations throw." << std::endl;

		// 3 parameters, each with correlations that form a non-PD matrix:
		// rho_12 = 0.9, rho_13 = 0.9, rho_23 = -0.9
		utility::tag::TagCOP perttag( tagptr_from_string(
			"<PerturbBundle name=\"pert\" default_perturbation_type=\"gaussian\">\n"
			"<Helix helix_index=\"1\" r0_perturbation=\"1.0\" omega0_perturbation=\"1.0\" delta_omega0_perturbation=\"1.0\" />\n"
			"<Correlation helix1=\"1\" param1=\"r0\" helix2=\"1\" param2=\"omega0\" correlation=\"0.9\" />\n"
			"<Correlation helix1=\"1\" param1=\"r0\" helix2=\"1\" param2=\"delta_omega0\" correlation=\"0.9\" />\n"
			"<Correlation helix1=\"1\" param1=\"omega0\" helix2=\"1\" param2=\"delta_omega0\" correlation=\"-0.9\" />\n"
			"</PerturbBundle>") );

		basic::datacache::DataMap dummy_data;
		protocols::helical_bundle::PerturbBundle pertbundle;
		TS_ASSERT_THROWS_ANYTHING( pertbundle.parse_my_tag( perttag, dummy_data ) );
	}

};
