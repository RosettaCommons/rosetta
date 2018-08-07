// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/helical_bundle/PerturbBundleTests.cxxtest.hh
/// @brief  Unit tests for the PerturbBundle mover.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>

// Project Headers
#include <protocols/helical_bundle/MakeBundle.hh>
#include <protocols/helical_bundle/PerturbBundle.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.hh>
#include <protocols/helical_bundle/parameters/BundleParametersSet.hh>
#include <protocols/helical_bundle/BundleParametrizationCalculator.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>

// Protocols Headers
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

static basic::Tracer TR("PerturbBundleTests");


class PerturbBundleTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}


	void test_perturb_r0_uniform() {
		TR << "Starting PerturbBundleTests::test_perturb_r0_uniform()." << std::endl;
		TR << "If this test should fail, it means that perturbation of float-valued parameters by the PerturbBundle mover is not working properly." << std::endl;

		utility::tag::TagCOP maketag( tagptr_from_string("<MakeBundle name=\"make\">\n<Helix r0=\"15\" />\n<Helix r0=\"25\" />\n</MakeBundle>") );
		utility::tag::TagCOP perttag( tagptr_from_string("<PerturbBundle name=\"pert\" default_perturbation_type=\"uniform\" >\n<Helix helix_index=\"2\" r0_perturbation=\"1.0\" />\n</PerturbBundle>") );
		core::Real smallest(500.0);
		core::Real largest(0.0);
		basic::datacache::DataMap dummy_data;
		protocols::filters::Filters_map dummy_filters;
		protocols::moves::Movers_map dummy_movers;

		for ( core::Size i(1); i<=30; ++i ) {
			core::pose::Pose testpose;

			protocols::helical_bundle::MakeBundle makebundle;
			makebundle.parse_my_tag(maketag, dummy_data, dummy_filters, dummy_movers, testpose );
			protocols::helical_bundle::PerturbBundle pertbundle;
			pertbundle.parse_my_tag(perttag, dummy_data, dummy_filters, dummy_movers, testpose );
			makebundle.apply(testpose);

			{
				protocols::helical_bundle::parameters::BundleParametersSetCOP paramset( utility::pointer::dynamic_pointer_cast< protocols::helical_bundle::parameters::BundleParametersSet const >( testpose.conformation().parameters_set(1) ) );
				TS_ASSERT( paramset != nullptr );
				protocols::helical_bundle::parameters::BundleParametersCOP params1( utility::pointer::dynamic_pointer_cast< protocols::helical_bundle::parameters::BundleParameters const >( paramset->parameters(1) ) );
				protocols::helical_bundle::parameters::BundleParametersCOP params2( utility::pointer::dynamic_pointer_cast< protocols::helical_bundle::parameters::BundleParameters const >( paramset->parameters(2) ) );
				TS_ASSERT( params1 != nullptr );
				TS_ASSERT( params2 != nullptr );
				core::conformation::parametric::RealValuedParameterCOP r0_1( utility::pointer::dynamic_pointer_cast< core::conformation::parametric::RealValuedParameter const >( params1->parameter_cop( protocols::helical_bundle::BPC_r0 ) ) );
				core::conformation::parametric::RealValuedParameterCOP r0_2( utility::pointer::dynamic_pointer_cast< core::conformation::parametric::RealValuedParameter const >( params2->parameter_cop( protocols::helical_bundle::BPC_r0 ) ) );
				TS_ASSERT_DELTA( r0_1->value(), 15.0, 1e-6 );
				TS_ASSERT_DELTA( r0_2->value(), 25.0, 1e-6 );
				TR << "Before:\tr0_1=" << r0_1->value() << "\tr0_2=" << r0_2->value() << std::endl;
			}
			pertbundle.apply(testpose);
			{
				protocols::helical_bundle::parameters::BundleParametersSetCOP paramset( utility::pointer::dynamic_pointer_cast< protocols::helical_bundle::parameters::BundleParametersSet const >( testpose.conformation().parameters_set(1) ) );
				TS_ASSERT( paramset != nullptr );
				protocols::helical_bundle::parameters::BundleParametersCOP params1( utility::pointer::dynamic_pointer_cast< protocols::helical_bundle::parameters::BundleParameters const >( paramset->parameters(1) ) );
				protocols::helical_bundle::parameters::BundleParametersCOP params2( utility::pointer::dynamic_pointer_cast< protocols::helical_bundle::parameters::BundleParameters const >( paramset->parameters(2) ) );
				TS_ASSERT( params1 != nullptr );
				TS_ASSERT( params2 != nullptr );
				core::conformation::parametric::RealValuedParameterCOP r0_1( utility::pointer::dynamic_pointer_cast< core::conformation::parametric::RealValuedParameter const >( params1->parameter_cop( protocols::helical_bundle::BPC_r0 ) ) );
				core::conformation::parametric::RealValuedParameterCOP r0_2( utility::pointer::dynamic_pointer_cast< core::conformation::parametric::RealValuedParameter const >( params2->parameter_cop( protocols::helical_bundle::BPC_r0 ) ) );
				TS_ASSERT_DELTA( r0_1->value(), 15.0, 1e-6 );
				TS_ASSERT_LESS_THAN( r0_2->value(), 26.0001 );
				TS_ASSERT_LESS_THAN( 23.9999, r0_2->value() );
				if ( r0_2->value() < smallest ) smallest = r0_2->value();
				if ( r0_2->value() > largest ) largest = r0_2->value();
				TR << "After:\tr0_1=" << r0_1->value() << "\tr0_2=" << r0_2->value() << std::endl;
			}
		}
		//After 30 tries, there should be some range of parameter values:
		TS_ASSERT_LESS_THAN( smallest, 24.95 );
		TS_ASSERT_LESS_THAN( 25.05, largest );
	}

};
