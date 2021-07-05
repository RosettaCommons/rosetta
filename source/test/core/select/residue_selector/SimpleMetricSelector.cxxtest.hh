// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/SimpleMetricSelector.cxxtest.hh
/// @brief test suite for core::select::residue_selector::SimpleMetricSelector
/// @author Brian Coventry


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>

// Package headers
#include <core/select/residue_selector/SimpleMetricSelector.hh>
#include <core/simple_metrics/per_residue_metrics/PerResidueEnergyMetric.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <protocols/rosetta_scripts/XmlObjects.hh>

// Utility headers
#include <utility/pointer/memory.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers

#include <core/init_util.hh> // AUTO IWYU For core_init

using namespace core;
using namespace core::select::residue_selector;

static basic::Tracer TR("core.select.residue_selector.SimpleMetricSelectorTests");


class SimpleMetricSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void do_test( core::pose::Pose const & pose, scoring::ScoreFunctionOP const & sfxn, Real lb, Real ub ) {
		if ( utility::isnan( lb ) && utility::isnan( ub) ) return;

		TR << "Test: " << lb << " " << ub << std::endl;

		core::simple_metrics::per_residue_metrics::PerResidueEnergyMetricOP metric =
			utility::pointer::make_shared<core::simple_metrics::per_residue_metrics::PerResidueEnergyMetric>();
		metric->set_scorefunction( sfxn );

		bool no_nan = ! ( utility::isnan( lb ) || utility::isnan( ub ) );

		SimpleMetricSelector sel_inside( metric, lb, ub );
		SimpleMetricSelector sel_outside( metric, lb, ub, no_nan );

		ResidueSubset sub_inside = sel_inside.apply( pose );
		ResidueSubset sub_outside = sel_outside.apply( pose );

		for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {

			Real value = pose.energies().residue_total_energy( seqpos );
			bool inside_lb = true;
			if ( ! utility::isnan( lb ) ) inside_lb = value >= lb;
			bool inside_ub = true;
			if ( ! utility::isnan( ub ) ) inside_ub = value <= ub;

			bool inside = inside_lb && inside_ub;

			TS_ASSERT( sub_inside[ seqpos ] == inside );

			if ( no_nan ) {
				TS_ASSERT( sub_outside[ seqpos ] != inside );
			}

		}

	}

	void test_selector() {

		core::pose::Pose pose = create_trpcage_ideal_pose();
		scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
		core::scoring::methods::EnergyMethodOptions emopts = sfxn->energy_method_options();
		emopts.hbond_options().decompose_bb_hb_into_pair_energies( true );
		sfxn->set_energy_method_options( emopts );

		sfxn->score( pose );

		utility::vector1<Real> choices { utility::get_undefined_real(), -10, -5, -2, -1, -0.5, 0, 0.5, 1, 2, 5, 10, utility::get_undefined_real() };

		for ( Size choice_lb = 1; choice_lb <= choices.size() - 1; choice_lb++ ) {
			for ( Size choice_ub = choice_lb + 1; choice_ub <= choices.size(); choice_ub++ ) {
				do_test( pose, sfxn, choices[ choice_lb ], choices[ choice_ub ] );
			}
		}
	}

	void test_xml() {


		protocols::rosetta_scripts::XmlObjectsCOP objs = protocols::rosetta_scripts::XmlObjects::create_from_string(
			"<SIMPLE_METRICS>\n"
			"<PerResidueEnergyMetric name=\"my_metric\"/>\n"
			"</SIMPLE_METRICS>\n"
			"<RESIDUE_SELECTORS>\n"
			"<SimpleMetricSelector name=\"sel1\" metric=\"my_metric\" lower_bound=\"\" upper_bound=\"15\" />\n"
			"<SimpleMetricSelector name=\"sel2\" lower_bound=\"15\" upper_bound=\"\" outside_bounds=\"true\" />\n"
			"</RESIDUE_SELECTORS>\n"
		);

		ResidueSelectorOP selector;

		selector = objs->get_residue_selector("sel1");
		SimpleMetricSelectorOP sel1 = std::dynamic_pointer_cast< SimpleMetricSelector >( selector );

		selector = objs->get_residue_selector("sel2");
		SimpleMetricSelectorOP sel2 = std::dynamic_pointer_cast< SimpleMetricSelector >( selector );

		TS_ASSERT( utility::isnan( sel1->lower_bound_ ) );
		TS_ASSERT_DELTA( sel1->upper_bound_, 15, 0.000001 );
		TS_ASSERT( ! sel1->outside_bounds_ );
		TS_ASSERT( sel1->metric_ );

		TS_ASSERT_DELTA( sel2->lower_bound_, 15, 0.000001 );
		TS_ASSERT( utility::isnan( sel2->upper_bound_ ) );
		TS_ASSERT( sel2->outside_bounds_ );
		TS_ASSERT( ! sel2->metric_ );


	}


};
