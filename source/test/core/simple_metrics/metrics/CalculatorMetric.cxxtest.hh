// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/metrics/CalculatorMetrics.cxxtest.hh
/// @brief  test for the Calculator Metric
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/rosettascripts.hh>
#include <core/simple_metrics/test_classes.hh>
#include <test/util/schema_utilities.hh>
#include <protocols/parser/SimpleMetricLoaderCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/extra_pose_info_util.hh>

#include <core/simple_metrics/metrics/CalculatorMetric.hh>
#include <core/simple_metrics/metrics/CalculatorMetricCreator.hh>

// Utility Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("core.simple_metrics.metrics.CalculatorMetric.cxxtest");

// --------------- Test Class --------------- //

class CalculatorMetricsTests : public CxxTest::TestSuite {

private:
	core::pose::PoseOP testpose_;
public:

	void setUp() {
		core_init();

		testpose_ = create_twores_1ubq_poseop(); // Identity doesn't matter
	}

	void tearDown() {
	}

	void test_calculatormetric() {
		using namespace core::simple_metrics;

		metrics::CalculatorMetric cf(" t1 = exp(-E1/kT); t2 = exp(-E2/kT); t1/( t1 + t2 ) " );

		cf.add_constant("kT", 0.6 );
		cf.add_simple_metric("E1", utility::pointer::make_shared< TestRealMetric >( -2 ) );
		cf.add_simple_metric("E2", utility::pointer::make_shared< TestRealMetric >( -1 ) );

		//default 0 threshold
		TS_ASSERT_DELTA( cf.calculate(*testpose_), 0.84113089511, 0.0001  );
	}

	void test_schema() {
		using namespace core::simple_metrics;

		utility::tag::XMLSchemaDefinition xsd;
		protocols::parser::SimpleMetricLoaderCreator().provide_xml_schema(xsd);

		metrics::CalculatorMetricCreator creator;
		check_if_tag_validates< metrics::CalculatorMetric >(
			// Taking advantage of C++11 raw string literals.
			R"xml(<CalculatorMetric name="test" equation="(min(a, b+b2, c)/(a*b*c-d)) + e + f" >
				<var name="a" metric="alpha" />
				<Var name="b" metric="beta" />
				<TestRealMetric name="b2" value="1" />
				<Var name="c" metric="delta" />
				<VAR name="d" value="-4.0" />
				<VAR name="e" reported="ten" />
				<TestRealMetric name="f" value="5.6" />
			</CalculatorMetric>
			)xml",
			xsd,
			creator.keyname(),
			complex_type_name_for_simple_metric( creator.keyname() )
		);
	}

	void test_parsing() {
		using namespace core::simple_metrics;

		basic::datacache::DataMap data;

		TestRealMetricOP sf1( new TestRealMetric( -1) );
		TestRealMetricOP sf2( new TestRealMetric( -2) );
		TestRealMetricOP sf3( new TestRealMetric( -3) );

		data["SimpleMetric"]["alpha"] = sf1;
		data["SimpleMetric"]["beta"] = sf2;
		data["SimpleMetric"]["delta"] = sf3;

		metrics::CalculatorMetric testmetric;
		TagCOP tag = tagptr_from_string(
			// Taking advantage of C++11 raw string literals.
			R"xml(<CalculatorMetric name="test" equation="(min(a, b+b2, c)/(a*b*c-d)) + e + f" >
				<var name="a" metric="alpha" />
				<Var name="b" metric="beta" />
				<TestRealMetric name="b2" value="1" />
				<Var name="c" metric="delta" />
				<VAR name="d" value="-4.0" />
				<VAR name="e" reported="ten" />
				<TestRealMetric name="f" value="5.6" />
			</CalculatorMetric>
			)xml"
		);

		testmetric.parse_my_tag( tag, data );
		setPoseExtraScore(*testpose_, "ten", 10.0);

		TS_ASSERT_EQUALS( testmetric.calculate( *testpose_), -3.0/(-1.0*-2.0*-3.0 - -4.0) + 10 + 5.6 );
	}

};
