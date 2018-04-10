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
#include <core/simple_metrics/StringMetric.hh>
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/IntegerMetric.hh>

#include <core/simple_metrics/CompositeStringMetric.hh>
#include <core/simple_metrics/CompositeRealMetric.hh>
#include <core/simple_metrics/CompositeIntegerMetric.hh>

#include <core/simple_metrics/test_classes.hh>

#include <core/simple_metrics/metrics/SelectedResiduesMetric.hh>
#include <core/simple_metrics/metrics/SelectedResiduesPyMOLMetric.hh>
#include <core/simple_metrics/metrics/TimingProfileMetric.hh>
#include <core/simple_metrics/metrics/SequenceMetric.hh>
#include <core/simple_metrics/metrics/SecondaryStructureMetric.hh>
#include <core/simple_metrics/metrics/SasaMetric.hh>

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
using namespace core::pose;
using namespace protocols::antibody::residue_selector;
using namespace protocols::antibody;
using namespace core::scoring;

class SimpleMetricTests : public CxxTest::TestSuite {
	//Define Variables

public:
	Pose pose; //Empty pose
	Pose ab_pose;

	void setUp(){
		core_init();

	}

	void test_string_metric() {
		TestStringMetric tester = TestStringMetric();
		TS_ASSERT( tester.calculate( pose ) == "TESTING");
		TS_ASSERT( tester.metric() == "SomeString");
		TS_ASSERT( tester.type() == "StringMetric");
		TS_ASSERT( tester.get_metric_names()[1] == "SomeString");

		tester.apply( pose );
		TS_ASSERT ( hasPoseExtraScore_str( pose, "SomeString" ) );

		std::string value;
		bool present = getPoseExtraScore( pose, "SomeString", value );
		TS_ASSERT( present );
		TS_ASSERT( value == "TESTING");

		tester.apply( pose, "prefix_", "_suffix");
		TS_ASSERT ( hasPoseExtraScore_str( pose, "prefix_SomeString_suffix" ) );
		present = getPoseExtraScore( pose, "prefix_SomeString_suffix", value );
		TS_ASSERT( present );
		TS_ASSERT( value == "TESTING");
	}

	void test_real_metric() {
		TestRealMetric tester = TestRealMetric();
		TS_ASSERT( tester.calculate( pose ) == 1.0);
		TS_ASSERT( tester.metric() == "SomeReal");
		TS_ASSERT( tester.type() == "RealMetric");
		TS_ASSERT( tester.get_metric_names()[1] == "SomeReal");

		tester.apply( pose );
		TS_ASSERT ( hasPoseExtraScore( pose, "SomeReal" ) );

		core::Real value;
		bool present = getPoseExtraScore( pose, "SomeReal", value );
		TS_ASSERT( present );

		TS_ASSERT( value == 1.0 );

		tester.apply( pose, "prefix_", "_suffix");
		TS_ASSERT ( hasPoseExtraScore( pose, "prefix_SomeReal_suffix" ) );
		present = getPoseExtraScore( pose, "prefix_SomeReal_suffix", value );
		TS_ASSERT( present );
		TS_ASSERT( value == 1.0 );

	}

	void test_integer_metric() {
		TestIntegerMetric tester = TestIntegerMetric();
		TS_ASSERT( tester.calculate( pose ) == 1);
		TS_ASSERT( tester.metric() == "SomeInteger");
		TS_ASSERT( tester.type() == "IntegerMetric");
		TS_ASSERT( tester.get_metric_names()[1] == "SomeInteger");

		tester.apply( pose );
		TS_ASSERT ( hasPoseExtraScore_int( pose, "SomeInteger" ) );

		int value;
		bool present = getPoseExtraScore( pose, "SomeInteger", value );
		TS_ASSERT( present );
		TS_ASSERT( value == 1 );

		tester.apply( pose, "prefix_", "_suffix");
		TS_ASSERT ( hasPoseExtraScore_int( pose, "prefix_SomeInteger_suffix" ) );
		present = getPoseExtraScore( pose, "prefix_SomeInteger_suffix", value );
		TS_ASSERT( present );
		TS_ASSERT( value == 1 );

	}

	void test_composite_string_metric() {
		TestCompositeStringMetric tester = TestCompositeStringMetric();
		std::map< std::string, std::string > values = tester.calculate( pose );
		utility::vector1< std::string > names = tester.get_metric_names();
		TS_ASSERT(names.size() == 2);
		TS_ASSERT( values["s_data1"] == "value1");
		TS_ASSERT( values["s_data2"] == "value2");

		TS_ASSERT( tester.metric() == "SomeCompositeString");
		TS_ASSERT( tester.type() == "CompositeStringMetric");

		tester.apply( pose );
		tester.apply( pose, "prefix_", "_suffix");

		for ( auto & name : names ) {

			TS_ASSERT ( hasPoseExtraScore_str( pose, name+"_SomeCompositeString" ) );

			std::string value;
			bool present = getPoseExtraScore( pose, name+"_SomeCompositeString", value );
			TS_ASSERT( present );
			TS_ASSERT( value ==  values[name] );


			TS_ASSERT ( hasPoseExtraScore_str( pose, "prefix_"+name+"_SomeCompositeString_suffix" ) );
			present = getPoseExtraScore( pose, "prefix_"+name+"_SomeCompositeString_suffix", value );
			TS_ASSERT( present );
			TS_ASSERT( value == values[name] );
		}

	}

	void test_composite_real_metric() {
		TestCompositeRealMetric tester = TestCompositeRealMetric();
		std::map< std::string, core::Real > values = tester.calculate( pose );
		utility::vector1< std::string > names = tester.get_metric_names();
		TS_ASSERT(names.size() == 2);
		TS_ASSERT( values["r_data1"] == 1.0);
		TS_ASSERT( values["r_data2"] == 2.0);

		TS_ASSERT( tester.metric() == "SomeCompositeReal");
		TS_ASSERT( tester.type() == "CompositeRealMetric");

		tester.apply( pose );
		tester.apply( pose, "prefix_", "_suffix");

		for ( auto & name : names ) {

			TS_ASSERT ( hasPoseExtraScore( pose, name+"_SomeCompositeReal" ) );

			core::Real value;
			bool present = getPoseExtraScore( pose, name+"_SomeCompositeReal", value );
			TS_ASSERT( present );
			TS_ASSERT( value ==  values[name] );


			TS_ASSERT ( hasPoseExtraScore( pose, "prefix_"+name+"_SomeCompositeReal_suffix" ) );
			present = getPoseExtraScore( pose, "prefix_"+name+"_SomeCompositeReal_suffix", value );
			TS_ASSERT( present );
			TS_ASSERT( value == values[name] );
		}

	}

	void test_composite_integer_metric() {
		TestCompositeIntegerMetric tester = TestCompositeIntegerMetric();
		std::map< std::string, int > values = tester.calculate( pose );
		utility::vector1< std::string > names = tester.get_metric_names();
		TS_ASSERT(names.size() == 2);
		TS_ASSERT( values["i_data1"] == 1);
		TS_ASSERT( values["i_data2"] == 2);

		TS_ASSERT( tester.metric() == "SomeCompositeInteger");
		TS_ASSERT( tester.type() == "CompositeIntegerMetric");

		tester.apply( pose );
		tester.apply( pose, "prefix_", "_suffix");

		for ( auto & name : names ) {

			TS_ASSERT ( hasPoseExtraScore_int( pose, name+"_SomeCompositeInteger" ) );

			int value;
			bool present = getPoseExtraScore( pose, name+"_SomeCompositeInteger", value );
			TS_ASSERT( present );
			TS_ASSERT( value ==  values[name] );


			TS_ASSERT ( hasPoseExtraScore_int( pose, "prefix_"+name+"_SomeCompositeInteger_suffix" ) );
			present = getPoseExtraScore( pose, "prefix_"+name+"_SomeCompositeInteger_suffix", value );
			TS_ASSERT( present );
			TS_ASSERT( value == values[name] );
		}

	}

	void test_utility_metrics() {
		core::import_pose::pose_from_file(ab_pose, "core/simple_metrics/2r0l_1_1.pdb", core::import_pose::PDB_file);
		AntibodyInfoOP ab_info =  AntibodyInfoOP( new AntibodyInfo(ab_pose, AHO_Scheme, North));
		CDRResidueSelectorOP cdr_selector = CDRResidueSelectorOP( new CDRResidueSelector( ab_info));
		cdr_selector->set_cdr( l1 );

		//SelectedResidues
		SelectedResiduesMetric selected_residues = SelectedResiduesMetric( cdr_selector );
		std::string pdb_nums = selected_residues.calculate( ab_pose );
		std::string pdb_nums_correct = "L24,L25,L26,L27,L28,L29,L38,L39,L40,L41,L42";
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
		AntibodyInfoOP ab_info =  AntibodyInfoOP( new AntibodyInfo(ab_pose, AHO_Scheme, North));
		CDRResidueSelectorOP cdr_selector = CDRResidueSelectorOP( new CDRResidueSelector( ab_info));
		cdr_selector->set_cdr( l1 );

		SecondaryStructureMetric ss_metric = SecondaryStructureMetric( cdr_selector );
		std::string l1_ss = ss_metric.calculate( ab_pose );
		std::string l1_correct = "EELLLLLLLEE";
		TS_ASSERT_EQUALS( l1_ss, l1_correct);
	}

	void test_seq_metric() {
		core::import_pose::pose_from_file(ab_pose, "core/simple_metrics/2r0l_1_1.pdb", core::import_pose::PDB_file);
		AntibodyInfoOP ab_info =  AntibodyInfoOP( new AntibodyInfo(ab_pose, AHO_Scheme, North));
		CDRResidueSelectorOP cdr_selector = CDRResidueSelectorOP( new CDRResidueSelector( ab_info));
		cdr_selector->set_cdr( l1 );

		SequenceMetric seq_metric = SequenceMetric( cdr_selector );
		std::string seq = seq_metric.calculate( ab_pose );
		std::string seq_correct = "DIQMTQSPSSL";
		TS_ASSERT_EQUALS( seq, seq_correct);
	}

	void test_sasa_metric() {
		core::import_pose::pose_from_file(ab_pose, "core/simple_metrics/2r0l_1_1.pdb", core::import_pose::PDB_file);
		AntibodyInfoOP ab_info =  AntibodyInfoOP( new AntibodyInfo(ab_pose, AHO_Scheme, North));
		CDRResidueSelectorOP cdr_selector = CDRResidueSelectorOP( new CDRResidueSelector( ab_info));
		cdr_selector->set_cdr( l1 );

		SasaMetric sasa_metric = SasaMetric ( );
		core::Real sas = sasa_metric.calculate( ab_pose );
		TS_ASSERT( sas > 496.5 );

		sasa_metric.set_residue_selector( cdr_selector );
		sas = sasa_metric.calculate( ab_pose );
		TS_ASSERT_DELTA( sas, 496.49, 1.0 );

	}


	void tearDown(){

	}







};
