// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/simple_metrics/SimpleMetricTests.cxxtest.hh
/// @brief  Test suite for core simple metrics.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/simple_metrics/util.hh>
#include <core/simple_metrics/SimpleMetricData.hh>
#include <core/simple_metrics/StringMetric.hh>
#include <core/simple_metrics/RealMetric.hh>

#include <core/simple_metrics/CompositeStringMetric.hh>
#include <core/simple_metrics/CompositeRealMetric.hh>

#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/PerResidueStringMetric.hh>

#include <core/simple_metrics/test_classes.hh>
#include <core/simple_metrics/test_classes.fwd.hh>

#include <core/simple_metrics/metrics/SelectedResiduesMetric.hh>
#include <core/simple_metrics/metrics/SelectedResiduesPyMOLMetric.hh>
#include <core/simple_metrics/metrics/TimingProfileMetric.hh>
#include <core/simple_metrics/metrics/SequenceMetric.hh>
#include <core/simple_metrics/metrics/SecondaryStructureMetric.hh>
#include <core/simple_metrics/metrics/SasaMetric.hh>
#include <core/simple_metrics/metrics/ResidueSummaryMetric.hh>
#include <core/simple_metrics/metrics/InteractionEnergyMetric.hh>
#include <core/simple_metrics/per_residue_metrics/PerResidueGlycanLayerMetric.hh>
#include <core/simple_metrics/per_residue_metrics/PerResidueClashMetric.hh>
#include <core/select/residue_selector/GlycanResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/io/raw_data/ScoreMap.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//
#include <protocols/antibody/residue_selector/CDRResidueSelector.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/AntibodyInfo.hh>

// Utility, etc Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <chrono>
#include <thread>

static basic::Tracer TR("SimpleMetricTests");

using namespace core::simple_metrics;
using namespace core::simple_metrics::metrics;
using namespace core::simple_metrics::per_residue_metrics;
using namespace core::pose;
using namespace protocols::antibody::residue_selector;
using namespace protocols::antibody;
using namespace core::scoring;
using namespace core::select;
using namespace core::select::residue_selector;
using namespace core::io::raw_data;

class SimpleMetricTests : public CxxTest::TestSuite {
	//Define Variables

public:
	Pose pose; //Empty pose
	Pose ab_pose;
	Pose glycan_pose;

	void setUp(){
		core_init();
		core_init_with_additional_options( "-include_sugars" );
	}

	void test_string_metric() {
		TestStringMetric tester = TestStringMetric();
		TS_ASSERT( tester.calculate( pose ) == "TESTING");
		TS_ASSERT( tester.metric() == "SomeString");
		TS_ASSERT( tester.simple_metric_type() == "StringMetric");
		TS_ASSERT( tester.get_metric_names()[1] == "SomeString");

		tester.apply( pose );

		TS_ASSERT ( has_sm_data( pose ));

		std::string value;
		bool present = get_sm_data(pose)->get_value( "SomeString", value );
		TS_ASSERT( present );
		TS_ASSERT( value == "TESTING");

		tester.apply( pose, "prefix_", "_suffix");
		present = get_sm_data( pose)->get_value("prefix_SomeString_suffix", value );
		TS_ASSERT( present );
		TS_ASSERT( value == "TESTING");

		//Test ScoreMap here
		std::map< std::string, std::string > scores = ScoreMap::get_arbitrary_string_data_from_pose( pose );
		TS_ASSERT( scores.count( "prefix_SomeString_suffix" ));
		TS_ASSERT( scores.at("prefix_SomeString_suffix") == "TESTING");

		//Test Cached Calculate function.
		TS_ASSERT_THROWS_NOTHING(tester.cached_calculate(pose, true, "prefix_", "_suffix", true /*fail_on_no_cache*/));
		value = tester.cached_calculate(pose, true, "prefix_", "_suffix", true);
		TS_ASSERT( value == "TESTING");

		//Test clearing data
		clear_sm_data( pose );
		TS_ASSERT( has_sm_data( pose )); //Class is still present - we do not kill it.
		present = get_sm_data( pose )->get_value("prefix_SomeString_suffix", value );
		TS_ASSERT( ! present );

	}

	void test_real_metric() {
		TestRealMetric tester = TestRealMetric();
		tester.set_custom_type("sometype");  //Test CustomType

		TS_ASSERT( tester.calculate( pose ) == 1.0);
		TS_ASSERT( tester.metric() == "SomeReal");
		TS_ASSERT( tester.simple_metric_type() == "RealMetric");
		TS_ASSERT( tester.get_metric_names()[1] == "SomeReal");

		tester.apply( pose );
		TS_ASSERT ( has_sm_data( pose ) );

		core::Real value;
		bool present = get_sm_data(pose)->get_value( "sometype_SomeReal", value );
		TS_ASSERT( present );

		TS_ASSERT( value == 1.0 );

		tester.apply( pose, "prefix_", "_suffix");
		present = get_sm_data( pose )->get_value("prefix_sometype_SomeReal_suffix", value );
		TS_ASSERT( present );
		TS_ASSERT( value == 1.0 );

		//Test ScoreMap here
		std::map< std::string, core::Real > scores = ScoreMap::get_arbitrary_score_data_from_pose( pose );
		TS_ASSERT( scores.count( "prefix_sometype_SomeReal_suffix" ) );
		TS_ASSERT( scores.at("prefix_sometype_SomeReal_suffix") == 1.0 );

		TS_ASSERT_THROWS_NOTHING(tester.cached_calculate(pose, true /*use_cache*/, "prefix_", "_suffix", true /*fail_on_no_cache*/));
		value = tester.cached_calculate(pose, true /*use_cache*/, "prefix_", "_suffix", true /*fail_on_no_cache*/);
		TS_ASSERT( value == 1.0 );

		//Test clearing all data
		clear_sm_data( pose );
		present = get_sm_data( pose )->get_value("prefix_sometype_SomeReal_suffix", value );
		TS_ASSERT( !present );
	}

	void test_composite_string_metric() {
		TestCompositeStringMetric tester = TestCompositeStringMetric();
		std::map< std::string, std::string > values = tester.calculate( pose );
		utility::vector1< std::string > names = tester.get_metric_names();
		TS_ASSERT(names.size() == 2);
		TS_ASSERT( values["s_data1"] == "value1");
		TS_ASSERT( values["s_data2"] == "value2");

		TS_ASSERT( tester.metric() == "SomeCompositeString");
		TS_ASSERT( tester.simple_metric_type() == "CompositeStringMetric");

		tester.apply( pose );
		tester.apply( pose, "prefix_", "_suffix");
		TS_ASSERT( has_sm_data(pose));

		bool present = get_sm_data(pose)->get_value("prefix_SomeCompositeString_suffix", values );
		TS_ASSERT( present );
		TS_ASSERT( values["s_data1"] == "value1");
		TS_ASSERT( values["s_data2"] == "value2");

		TS_ASSERT_THROWS_NOTHING(tester.cached_calculate(pose, true /*use_cache*/, "prefix_", "_suffix", true /*fail_on_no_cache*/));
		values = tester.cached_calculate(pose, true /*use_cache*/, "prefix_", "_suffix", true /*fail_on_no_cache*/);
		TS_ASSERT( values["s_data1"] == "value1");
		TS_ASSERT( values["s_data2"] == "value2");

		//Test clearing of all data
		clear_sm_data( pose );
		present = get_sm_data(pose)->get_value("prefix_SomeCompositeString_suffix", values );
		TS_ASSERT( ! present );

	}

	void test_composite_real_metric() {
		TestCompositeRealMetric tester = TestCompositeRealMetric();
		std::map< std::string, core::Real > values = tester.calculate( pose );
		utility::vector1< std::string > names = tester.get_metric_names();
		TS_ASSERT(names.size() == 2);
		TS_ASSERT( values["r_data1"] == 1.0);
		TS_ASSERT( values["r_data2"] == 2.0);

		TS_ASSERT( tester.metric() == "SomeCompositeReal");
		TS_ASSERT( tester.simple_metric_type() == "CompositeRealMetric");

		tester.apply( pose );
		tester.apply( pose, "prefix_", "_suffix");
		TS_ASSERT(has_sm_data(pose));

		bool present = get_sm_data( pose )->get_value( "prefix_SomeCompositeReal_suffix", values );
		TS_ASSERT( present );
		TS_ASSERT( values["r_data1"] == 1.0);
		TS_ASSERT( values["r_data2"] == 2.0);

		TS_ASSERT_THROWS_NOTHING(tester.cached_calculate(pose, true /*use_cache*/, "prefix_", "_suffix", true /*fail_on_no_cache*/));
		values = tester.cached_calculate(pose, true /*use_cache*/, "prefix_", "_suffix", true /*fail_on_no_cache*/);
		TS_ASSERT( values["r_data1"] == 1.0);
		TS_ASSERT( values["r_data2"] == 2.0);

		clear_sm_data( pose );
		present = get_sm_data( pose )->get_value( "prefix_SomeCompositeReal_suffix", values );
		TS_ASSERT( ! present );

	}

	void test_per_residue_real_metric() {
		using namespace core::select::residue_selector;
		GlycanResidueSelectorOP selector = utility::pointer::make_shared< GlycanResidueSelector >();

		TestPerResidueRealMetric tester = TestPerResidueRealMetric();
		tester.set_residue_selector(selector); //Test setting the residue selector.
		std::map< core::Size, core::Real > values = tester.calculate( pose );
		utility::vector1< std::string > const names = tester.get_metric_names();
		TS_ASSERT(names.size() == 1);
		TS_ASSERT( values.at(1) == 1.0);
		TS_ASSERT( values.at(2) == 2.0);

		TS_ASSERT( tester.metric() == "SomePerResidueReal");
		TS_ASSERT( tester.simple_metric_type() == "PerResidueRealMetric");

		tester.apply( pose );
		tester.apply( pose, "prefix_", "_suffix");
		TS_ASSERT( has_sm_data(pose));

		bool present = get_sm_data(pose)->get_value( "prefix_SomePerResidueReal_suffix", values);
		TS_ASSERT(present);
		TS_ASSERT( values.at(1) == 1.0);
		TS_ASSERT( values.at(2) == 2.0);

		TS_ASSERT_THROWS_NOTHING(tester.cached_calculate(pose, true /*use_cache*/, "prefix_", "_suffix", true /*fail on missing cache*/, false /* use ref pose */));

		values = tester.cached_calculate(pose, true /*use_cache*/, "prefix_", "_suffix", true /*fail_on_no_cache*/, false /*refpose*/);
		TS_ASSERT( values.at(1) == 1.0);
		TS_ASSERT( values.at(2) == 2.0);

		clear_sm_data( pose );
		present = get_sm_data(pose)->get_value( "prefix_SomePerResidueReal_suffix", values);
		TS_ASSERT( ! present );

	}

	void test_per_residue_string_metric() {
		using namespace core::select::residue_selector;
		GlycanResidueSelectorOP selector = utility::pointer::make_shared< GlycanResidueSelector >();

		TestPerResidueStringMetric tester = TestPerResidueStringMetric();
		tester.set_residue_selector(selector); //Test setting the residue selector.
		std::map< core::Size, std::string > values = tester.calculate( pose );
		utility::vector1< std::string > const names = tester.get_metric_names();
		TS_ASSERT(names.size() == 1);
		TS_ASSERT( values.at(1) == "value1");
		TS_ASSERT( values.at(2) == "value2");

		TS_ASSERT( tester.metric() == "SomePerResidueString");
		TS_ASSERT( tester.simple_metric_type() == "PerResidueStringMetric");

		bool use_ref_pose = false;

		tester.apply( pose );
		tester.apply( pose, "prefix_", "_suffix");
		TS_ASSERT( has_sm_data(pose));

		bool present = get_sm_data(pose)->get_value( "prefix_SomePerResidueString_suffix", values);
		TS_ASSERT(present);
		TS_ASSERT( values.at(1) == "value1");
		TS_ASSERT( values.at(2) == "value2");

		TS_ASSERT_THROWS_NOTHING(tester.cached_calculate(pose, true /*use_cache*/, "prefix_", "_suffix", true /*fail_on_no_cache*/, use_ref_pose ));
		values = tester.cached_calculate(pose, true /*use_cache*/, "prefix_", "_suffix", true, use_ref_pose);
		TS_ASSERT( values.at(1) == "value1");
		TS_ASSERT( values.at(2) == "value2");

		clear_sm_data( pose );
		present = get_sm_data(pose)->get_value( "prefix_SomePerResidueString_suffix", values);
		TS_ASSERT( ! present );


	}

	void test_utility_metrics() {
		core::import_pose::pose_from_file(ab_pose, "core/simple_metrics/2r0l_1_1.pdb", core::import_pose::PDB_file);
		AntibodyInfoOP ab_info =  utility::pointer::make_shared< AntibodyInfo >(ab_pose, AHO_Scheme, North);
		CDRResidueSelectorOP cdr_selector = utility::pointer::make_shared< CDRResidueSelector >( ab_info);
		cdr_selector->set_cdr( l1 );

		//SelectedResidues
		SelectedResiduesMetric selected_residues = SelectedResiduesMetric( cdr_selector );
		std::string pdb_nums = selected_residues.calculate( ab_pose );
		std::string pdb_nums_correct = "24L,25L,26L,27L,28L,29L,38L,39L,40L,41L,42L";
		TS_ASSERT_EQUALS( pdb_nums, pdb_nums_correct);

		selected_residues.set_output_in_rosetta_num( true );
		std::string rosetta_nums = selected_residues.calculate( ab_pose );
		std::string rosetta_nums_correct = "24,25,26,27,28,29,30,31,32,33,34";
		TS_ASSERT_EQUALS( rosetta_nums, rosetta_nums_correct );

		//SelectedResiduesPyMOL
		SelectedResiduesPyMOLMetric selected_pymol = SelectedResiduesPyMOLMetric( cdr_selector );
		std::string pymol_sele = selected_pymol.calculate( ab_pose );
		std::string pymol_sele_correct = "select rosetta_sele, (chain L and resid 24,25,26,27,28,29,38,39,40,41,42)";
		TS_ASSERT_EQUALS( pymol_sele, pymol_sele_correct);

		//TimingProfileMetric
		TimingProfileMetric timing = TimingProfileMetric();
		std::this_thread::sleep_for (std::chrono::seconds(10));
		core::Real time_since_construction = timing.calculate( ab_pose );
		TS_ASSERT( time_since_construction > 0); //Variable time - hard to test exactly, - just make sure it is not zero.

	}

	void test_ss_metric() {
		core::import_pose::pose_from_file(ab_pose, "core/simple_metrics/2r0l_1_1.pdb", core::import_pose::PDB_file);
		AntibodyInfoOP ab_info =  utility::pointer::make_shared< AntibodyInfo >(ab_pose, AHO_Scheme, North);
		CDRResidueSelectorOP cdr_selector = utility::pointer::make_shared< CDRResidueSelector >( ab_info);
		cdr_selector->set_cdr( l1 );

		SecondaryStructureMetric ss_metric = SecondaryStructureMetric( cdr_selector );
		std::string l1_ss = ss_metric.calculate( ab_pose );
		std::string l1_correct = "EELLLLLLLEE";
		TS_ASSERT_EQUALS( l1_ss, l1_correct);
	}

	void test_seq_metric() {
		core::import_pose::pose_from_file(ab_pose, "core/simple_metrics/2r0l_1_1.pdb", core::import_pose::PDB_file);
		AntibodyInfoOP ab_info =  utility::pointer::make_shared< AntibodyInfo >(ab_pose, AHO_Scheme, North);
		CDRResidueSelectorOP cdr_selector = utility::pointer::make_shared< CDRResidueSelector >( ab_info);
		cdr_selector->set_cdr( l1 );

		SequenceMetric seq_metric = SequenceMetric( cdr_selector );
		std::string seq = seq_metric.calculate( ab_pose );
		std::string seq_correct = "RASQDVSTAVA";
		TS_ASSERT_EQUALS( seq, seq_correct);
	}

	void test_sasa_metric() {
		core::import_pose::pose_from_file(ab_pose, "core/simple_metrics/2r0l_1_1.pdb", core::import_pose::PDB_file);
		AntibodyInfoOP ab_info =  utility::pointer::make_shared< AntibodyInfo >(ab_pose, AHO_Scheme, North);
		CDRResidueSelectorOP cdr_selector = utility::pointer::make_shared< CDRResidueSelector >( ab_info);
		cdr_selector->set_cdr( l1 );

		SasaMetric sasa_metric = SasaMetric ( );
		core::Real const sas = sasa_metric.calculate( ab_pose );
		TR << "Total sasa: " << sas << std::endl;
		TS_ASSERT( sas > 496.5 );

		sasa_metric.set_residue_selector( cdr_selector );
		core::Real const sas_subset = sasa_metric.calculate( ab_pose );
		TR << "Subset sasa: " << sas_subset << std::endl;
		TS_ASSERT_DELTA( sas_subset, 496.49, 1.0 );

		sasa_metric.set_sasa_metric_mode( core::scoring::sasa::SasaMethodHPMode::POLAR_SASA );
		core::Real const sas_subset_polar( sasa_metric.calculate( ab_pose ) );
		TR << "Subset polar sasa: " << sas_subset_polar << std::endl;
		TS_ASSERT_LESS_THAN( sas_subset_polar, sas_subset );
		TS_ASSERT_DELTA( sas_subset_polar, 270.297, 1.0 );

		sasa_metric.set_sasa_metric_mode( core::scoring::sasa::SasaMethodHPMode::HYDROPHOBIC_SASA );
		core::Real const sas_subset_apolar( sasa_metric.calculate( ab_pose ) );
		TR << "Subset apolar sasa: " << sas_subset_apolar << std::endl;
		TS_ASSERT_LESS_THAN( sas_subset_apolar, sas_subset );
		TS_ASSERT_DELTA( sas_subset_apolar, 226.196, 1.0 );

		TS_ASSERT_DELTA( sas_subset_apolar + sas_subset_polar, sas_subset, 0.01 );
	}

	void test_glycan_layer_metric() {
		pose_from_file( glycan_pose, "core/chemical/carbohydrates/gp120_2glycans_man5.pdb" , core::import_pose::PDB_file);
		PerResidueGlycanLayerMetric layer_metric = PerResidueGlycanLayerMetric();
		std::map< core::Size, core::Real > result = layer_metric.calculate(glycan_pose);

		for ( auto & pair : result ) {
			core::Size layer = static_cast<core::Size>(pair.second);
			TS_ASSERT(layer == glycan_pose.glycan_tree_set()->get_distance_to_start(pair.first));
		}
	}

	void test_residue_summary_metric() {
		using namespace core::select::residue_selector;


		TestPerResidueRealMetricOP tester = TestPerResidueRealMetricOP( new TestPerResidueRealMetric());
		ResidueSummaryMetric summary_metric = ResidueSummaryMetric(tester);

		//std::map< core::Size, core::Real > values = tester.calculate( pose );
		//TS_ASSERT(names.size() == 1);
		//TS_ASSERT( values.at(1) == 1.0);
		//TS_ASSERT( values.at(2) == 2.0);


		summary_metric.set_action( sum );
		TS_ASSERT_DELTA( summary_metric.calculate(pose), 3.0, .01);

		summary_metric.set_action( mean );
		TS_ASSERT_DELTA( summary_metric.calculate(pose), 1.5, .01);

		summary_metric.set_action( n_res_eq);
		summary_metric.set_action_value( 1.0 );
		TS_ASSERT( summary_metric.calculate(pose) == 1.0 );

		summary_metric.set_action( n_res_ne);
		TS_ASSERT( summary_metric.calculate(pose) == 1.0 );

		summary_metric.set_action( n_res_gt );
		summary_metric.set_action_value( 0.0 );
		TS_ASSERT( summary_metric.calculate(pose) == 2.0 );

		summary_metric.set_action( n_res_gt_or_eq );
		summary_metric.set_action_value( 1.0 );
		TS_ASSERT( summary_metric.calculate(pose ) == 2.0 );

		summary_metric.set_action( n_res_lt_or_eq );
		summary_metric.set_action_value( 2.0 );
		TS_ASSERT( summary_metric.calculate(pose ) == 2.0);

		summary_metric.set_action( n_res_lt );
		TS_ASSERT( summary_metric.calculate(pose) == 1.0 );


		/// Test Caching
		tester->apply(ab_pose);
		summary_metric.set_use_cached_data(true);
		summary_metric.set_fail_on_missing_cache(true);

		TS_ASSERT( summary_metric.calculate(ab_pose) == 1.0 );


	}

	void test_interaction_energy_metric(){
		core::import_pose::pose_from_file(ab_pose, "core/simple_metrics/2r0l_1_1.pdb", core::import_pose::PDB_file);
		ScoreFunctionOP default_scorefxn = get_score_function();
		default_scorefxn->score(ab_pose);
		AntibodyInfoOP ab_info =  AntibodyInfoOP( new AntibodyInfo(ab_pose, AHO_Scheme, North));
		CDRResidueSelectorOP cdr_selector = CDRResidueSelectorOP( new CDRResidueSelector( ab_info));
		cdr_selector->set_cdr( l3 );

		NotResidueSelectorOP not_L3 = NotResidueSelectorOP( new NotResidueSelector( cdr_selector));

		InteractionEnergyMetric metric = InteractionEnergyMetric( cdr_selector, not_L3 );
		TS_ASSERT( metric.calculate(ab_pose) < 0);  //In case of score function changes, we just check this is indeed NOT zero


		ScoreFunctionOP scorefxn = ScoreFunctionOP( new ScoreFunction());
		scorefxn->set_weight(fa_rep, 1.0);
		scorefxn->score(ab_pose);

		utility::vector1< ScoreType > types;
		types.push_back(fa_rep);

		metric.set_include_only_scoretypes( types );
		TS_ASSERT( metric.calculate(ab_pose) >= 12.0);


	}

	void test_per_residue_clash_metric(){
		core::import_pose::pose_from_file(glycan_pose, "core/simple_metrics/two_glycans.pdb", core::import_pose::PDB_file);
		GlycanResidueSelectorOP glycan_selector = GlycanResidueSelectorOP( new GlycanResidueSelector(148, true /*include_root*/));

		ScoreFunctionOP default_scorefxn = get_score_function();
		default_scorefxn->score(glycan_pose);

		NotResidueSelectorOP not_glycans = NotResidueSelectorOP( new NotResidueSelector( glycan_selector));

		PerResidueClashMetricOP metric = PerResidueClashMetricOP(new PerResidueClashMetric(glycan_selector, not_glycans));
		metric->set_use_soft_clash(true);
		metric->set_use_hydrogens(false);
		metric->set_soft_dampening(.33);

		ResidueSummaryMetric summary_metric = ResidueSummaryMetric( metric );
		summary_metric.set_action(n_res_gt);
		summary_metric.set_action_value(0.0);

		//No clashes
		TS_ASSERT(summary_metric.calculate(glycan_pose) == 0);

		metric->set_secondary_residue_selector(glycan_selector);
		summary_metric.set_metric(metric); //Not needed - but just to make sure.
		TS_ASSERT_DELTA(summary_metric.calculate(glycan_pose), 7.0, .01);


	}

	void tearDown(){

	}







};
