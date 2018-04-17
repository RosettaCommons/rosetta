// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/analysis/simple_metrics/SimpleMetricProtocolsTests.cxxtest.hh
/// @brief  Test suite for simple metrics in protocols.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/simple_metrics/test_classes.hh>
#include <protocols/simple_filters/SimpleMetricFilter.hh>
#include <protocols/analysis/simple_metrics/RMSDMetric.hh>
#include <protocols/analysis/simple_metrics/DihedralDistanceMetric.hh>
#include <protocols/analysis/simple_metrics/TotalEnergyMetric.hh>
#include <protocols/analysis/simple_metrics/CompositeEnergyMetric.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>

using namespace core::simple_metrics;
using namespace protocols::simple_filters;
using namespace protocols::analysis::simple_metrics;
using namespace core::pose;
using namespace core::select::residue_selector;
using namespace core::scoring;

static basic::Tracer TR("SimpleMetricProtocolsTests");

class SimpleMetricProtocolsTests : public CxxTest::TestSuite {
	//Define Variables

public:

	SimpleMetricFilter filter;
	//protocols::antibody::AntibodyInfoOP ab_info;
	Pose pose;

	void setUp(){
		core_init();

		core::import_pose::pose_from_file(pose, "protocols/antibody/2r0l_1_1.pdb", core::import_pose::PDB_file);
		//ab_info =  protocols::antibody::AntibodyInfoOP( new protocols::antibody::AntibodyInfo(ab_pose, AHO_Scheme, North));
		//CDRResidueSelectorOP cdr_selector = CDRResidueSelectorOP( new CDRResidueSelector( ab_info));


	}

	void tearDown(){

	}

	void test_string_metric(){
		TestStringMetricOP tester = TestStringMetricOP( new TestStringMetric());
		filter.set_simple_metric( tester );
		filter.set_match_string( "TESTING");
		filter.set_comparison_type( eq );

		TS_ASSERT( filter.apply(pose));

		filter.set_comparison_type( ne );
		TS_ASSERT( ! filter.apply( pose ));

	}

	void test_real_metric() {
		TestRealMetricOP tester = TestRealMetricOP( new TestRealMetric());
		filter.set_simple_metric( tester );

		filter.set_cutoff(0);

		filter.set_comparison_type( eq );
		TS_ASSERT( ! filter.apply( pose ));

		filter.set_comparison_type( ne );
		TS_ASSERT( filter.apply( pose ));

		filter.set_comparison_type( lt );
		TS_ASSERT( ! filter.apply( pose ));

		filter.set_comparison_type( gt );
		TS_ASSERT( filter.apply( pose ));

		filter.set_comparison_type( lt_or_eq );
		TS_ASSERT( ! filter.apply( pose ));

		filter.set_comparison_type( gt_or_eq );
		TS_ASSERT( filter.apply( pose ));

		filter.set_cutoff(1.0);
		filter.set_comparison_type( eq );
		TS_ASSERT( filter.apply( pose ));

		filter.set_comparison_type( lt_or_eq);
		TS_ASSERT( filter.apply( pose ));

		filter.set_comparison_type( gt_or_eq );
		TS_ASSERT( filter.apply( pose ));

	}

	void test_composite_string_metric() {
		TestCompositeStringMetricOP tester = TestCompositeStringMetricOP( new TestCompositeStringMetric());
		filter.set_simple_metric( tester );
		filter.set_comparison_type(eq);

		filter.set_composite_action("any");
		filter.set_match_string( "TESTING");
		TS_ASSERT( ! filter.apply( pose ));

		filter.set_match_string( "value1");
		TS_ASSERT( filter.apply( pose ));

		filter.set_match_string( "value2");
		TS_ASSERT( filter.apply( pose ) );

		filter.set_composite_action("all");
		TS_ASSERT( ! filter.apply( pose ));

		filter.set_composite_action("s_data1");
		TS_ASSERT( ! filter.apply( pose ));

		filter.set_composite_action("s_data2");
		TS_ASSERT( filter.apply( pose ));

		filter.set_comparison_type(ne);
		TS_ASSERT( ! filter.apply( pose ))

	}

	void test_composite_real_metric() {

		TestCompositeRealMetricOP tester = TestCompositeRealMetricOP( new TestCompositeRealMetric());
		filter.set_simple_metric( tester );
		filter.set_comparison_type( eq );

		filter.set_cutoff(0);
		filter.set_composite_action("any");

		TS_ASSERT( ! filter.apply( pose ));

		filter.set_comparison_type( lt );
		TS_ASSERT( ! filter.apply( pose ))

			filter.set_comparison_type( gt );
		TS_ASSERT( filter.apply( pose ));

		filter.set_composite_action("all");
		TS_ASSERT( filter.apply( pose ));

		filter.set_composite_action("r_data1");
		TS_ASSERT (filter.apply( pose ));

	}

	void test_total_energy_metric(){
		ScoreFunctionOP scorefxn = get_score_function();
		TotalEnergyMetric e_metric = TotalEnergyMetric();
		e_metric.set_scorefunction(scorefxn);

		TS_ASSERT_DELTA( e_metric.calculate( pose ), scorefxn->score( pose ), 1);
	}

	void test_composite_total_energy_metric() {
		ScoreFunctionOP scorefxn = get_score_function();
		CompositeEnergyMetric e_metric = CompositeEnergyMetric();
		e_metric.set_scorefunction(scorefxn);

		std::map< std::string, core::Real > energies = e_metric.calculate( pose );
		scorefxn->score( pose );

		EnergyMap weights( scorefxn->weights() );
		EnergyMap all_energies = ( weights * pose.energies().total_energies() );

		TS_ASSERT( energies.size() != 1 );

		//Test to make sure all non-zero weights are present
		for ( core::Size ii = 1; ii <= n_score_types; ++ii ) {
			ScoreType score_type = static_cast<ScoreType>(ii);

			if ( scorefxn->has_nonzero_weight( score_type ) ) {
				std::string keyname= name_from_score_type( score_type  );
				TS_ASSERT( energies.count( keyname) != 0 );

				TS_ASSERT_DELTA( energies[ keyname ], all_energies[ score_type ], .1);

			}
		}

	}

	void test_rmsd_metric() {
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueIndexSelectorOP selector1 = ResidueIndexSelectorOP( new ResidueIndexSelector( "1-10" ));
		ResidueIndexSelectorOP selector2 = ResidueIndexSelectorOP( new ResidueIndexSelector( "11-20" ));

		RMSDMetric rmsd_metric = RMSDMetric(trpcage);

		TS_ASSERT_EQUALS( rmsd_metric.calculate(trpcage), 0.0 );

		rmsd_metric.set_residue_selector( selector1 );
		rmsd_metric.set_residue_selector_reference( selector2 );
		rmsd_metric.set_rmsd_type( rmsd_protein_bb_ca );

		TS_ASSERT_DELTA( rmsd_metric.calculate( trpcage ), 13.0370, 0.01 );
	}

	void test_dihedral_metric() {
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueIndexSelectorOP selector1 = ResidueIndexSelectorOP( new ResidueIndexSelector( "1-10" ));
		ResidueIndexSelectorOP selector2 = ResidueIndexSelectorOP( new ResidueIndexSelector( "11-20" ));

		DihedralDistanceMetric dih_metric = DihedralDistanceMetric();
		dih_metric.set_comparison_pose( trpcage );

		TS_ASSERT_DELTA( dih_metric.calculate( trpcage ), 0.0, 0.1);

		dih_metric.set_residue_selector( selector1 );
		dih_metric.set_residue_selector_reference( selector2 );

		TS_ASSERT_DELTA( dih_metric.calculate( trpcage ), 76.838, 1.0);

	}

};
