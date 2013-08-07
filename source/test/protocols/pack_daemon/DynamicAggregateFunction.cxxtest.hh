// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pack_daemon/PackDaemon.cxxtest.hh
/// @brief  test suite for PackDaemon
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <protocols/pack_daemon/DynamicAggregateFunction.hh>
#include <protocols/pack_daemon/PackDaemon.hh>
// AUTO-REMOVED #include <protocols/pack_daemon/EntityCorrespondence.hh>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Core headers
// AUTO-REMOVED #include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerSet.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerSets.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Protocols headers
#include <protocols/genetic_algorithm/Entity.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <sstream>

//Auto Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>


//static basic::Tracer TR("PackDaemonTest.cxxtest");

using namespace core;
using namespace protocols::pack_daemon;
using namespace protocols::genetic_algorithm;
using namespace core::chemical;
using namespace core::scoring;
using namespace core::pack::rotamer_set;

class TestCalculator : public NPDPropCalculator
{
public:

	virtual
	core::Real
	calculate( core::pose::Pose const & /* p */ ) { return 1.0; };

};

class TestCalculatorCreator : public NPDPropCalculatorCreator
{
	virtual
	std::string
	calculator_name() const {return "TESTCALC"; }

	virtual
	NPDPropCalculatorOP
	new_calculator() const { return new TestCalculator; }
};


class DynamicAggregateFunctionTests : public CxxTest::TestSuite
{
public:
	typedef core::pose::PoseOP PoseOP;
	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef core::pack::task::TaskFactory TaskFactory;
	typedef protocols::pack_daemon::PackDaemon PackDaemon;

public:
	void setUp() {
		core_init();
	}

	void test_initialize_empty_DAF() {
		DaemonSetOP ds = new DaemonSet;
		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;

		std::istringstream iss( "FITNESS 3.0 + 4.0" );

		daf->initialize_from_input_file( ds, iss );

	}

	void test_daf_STATE_command() {

		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );

		std::istringstream iss( "STATE trpcage 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\nFITNESS trpcage + 2.0\n" );

		daf->initialize_from_input_file( ds, iss );

	}

	void test_daf_STATE_command_missing_secondary_resfile() {

		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );

		std::istringstream iss( "STATE trpcage 1l2y.pdb 1l2y.correspondence.txt\nFITNESS trpcage + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string expected_error_message = "Expected to read secondary resfile name in the DynamicAggregateFunction input file after reading STATE correspondence file on line 1\nSTATE trpcage 1l2y.pdb 1l2y.correspondence.txt";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}

	}

	void test_daf_STATE_command_missing_correspondence_file() {

		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );

		std::istringstream iss( "STATE trpcage 1l2y.pdb\nFITNESS trpcage + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string expected_error_message = "Expected to read secondary resfile name in the DynamicAggregateFunction input file after reading STATE correspondence file on line 1\nSTATE trpcage 1l2y.pdb";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}

	}

	void test_daf_STATE_command_missing_pdb_file() {

		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );

		std::istringstream iss( "STATE trpcage\nFITNESS trpcage + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << e.msg();
			std::string expected_error_message =
				"Expected to read pdb file name in the DynamicAggregateFunction input file after reading STATE pdb file on line 1\n"
				"STATE trpcage";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}

	}

	void test_daf_STATE_command_missing_secondary_resfile2() {

		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );

		std::istringstream iss( "STATE trpcage 1l2y.pdb 1l2y.correspondence.txt" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << e.msg();
			std::string expected_error_message = "Expected to read secondary resfile name in the DynamicAggregateFunction input file after reading STATE correspondence file on line 1\nSTATE trpcage 1l2y.pdb 1l2y.correspondence.txt";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}

	}

	void test_daf_STATE_command_missing_correspondence_file2() {

		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );

		std::istringstream iss( "STATE trpcage 1l2y.pdb" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << e.msg();
			std::string expected_error_message = "Expected to read secondary resfile name in the DynamicAggregateFunction input file after reading STATE correspondence file on line 1\nSTATE trpcage 1l2y.pdb";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}

	}

	void test_daf_STATE_command_missing_pdb_file2() {

		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );

		std::istringstream iss( "STATE trpcage" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << e.msg();
			std::string expected_error_message =
				"Expected to read pdb file name in the DynamicAggregateFunction input file after reading STATE pdb file on line 1\n"
				"STATE trpcage";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}

	}

	void test_daf_STATE_command_varname_illegal_name() {

		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );

		std::istringstream iss( "STATE max 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\nFITNESS trpcage + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << "test_daf_STATE_command_varname_bad" << std::endl;
			//std::cout << e.msg() << std::endl;
			std::string expected_error_message = "Illegal name for variable, 'max' in the STATE command on line 1\nSTATE max 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}

	}

	void test_daf_STATE_command_varname_function_name() {

		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );

		std::istringstream iss( "STATE pow 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\nFITNESS trpcage + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << "test_daf_STATE_command_varname_function_name" << std::endl;
			//std::cout << e.msg() << std::endl;
			//std::string expected_error_message = "Declaration of variable 'pow' in STATE conflicts with a function name.  Line 1\nSTATE pow 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile";
			std::string expected_error_message =
				"Declaration of variable 'pow' in STATE command conflicts with a function name.  Line 1\n"
				"STATE pow 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}

	}

	void test_daf_STATE_command_varname_duplicated() {

		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );

		std::istringstream iss( "STATE trpcage 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\nSTATE trpcage 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\nFITNESS trpcage + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << "test_daf_STATE_command_varname_duplicated" << std::endl;
			//std::cout << e.msg() << std::endl << std::endl;
			//std::string expected_error_message = "Variable name trpcage appears multiple times in the DynamicAggregateFunction file.\nFirst occurrance was found on line 1.  Second occurrance found while reading STATE command\nSTATE trpcage 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\nLine # 2";
			std::string expected_error_message =
				"Variable name trpcage appears multiple times in the DynamicAggregateFunction file.\n"
				"Second occurrance found while reading a STATE command\n"
				"STATE trpcage 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n"
				"Line # 2";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}

	}


	void test_daf_STATE_VECTOR_command() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );

		std::istringstream iss( "STATE_VECTOR trpcage trpcage.list\nFITNESS vmin( trpcage ) + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( true );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cout << e.msg();
			TS_ASSERT( false );
		}
	}

	void test_daf_STATE_VECTOR_command_missing_listfile() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );

		std::istringstream iss( "STATE_VECTOR trpcage\nFITNESS vmin( trpcage ) + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << e.msg();
			std::string expected_error_message = "Expected to read state-vector file name in the DynamicAggregateFunction input file after reading STATE_VECTOR variable name on line 1\nSTATE_VECTOR trpcage";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}
	}

	void test_daf_STATE_VECTOR_command_missing_varname() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );

		std::istringstream iss( "STATE_VECTOR\nFITNESS vmin( trpcage ) + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << e.msg();
			std::string expected_error_message = "Expected to read state-vector variable name in the DynamicAggregateFunction input file after reading STATE_VECTOR on line 1\nSTATE_VECTOR";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}
	}

	void test_daf_STATE_VECTOR_command_missing_listfile2() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );

		std::istringstream iss( "STATE_VECTOR trpcage" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << e.msg();
			std::string expected_error_message = "Expected to read state-vector file name in the DynamicAggregateFunction input file after reading STATE_VECTOR variable name on line 1\nSTATE_VECTOR trpcage";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}
	}

	void test_daf_STATE_VECTOR_command_missing_varname2() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );

		std::istringstream iss( "STATE_VECTOR" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << e.msg();
			std::string expected_error_message = "Expected to read state-vector variable name in the DynamicAggregateFunction input file after reading STATE_VECTOR on line 1\nSTATE_VECTOR";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}
	}

	void test_daf_STATE_VECTOR_command_illegal_varname() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );

		std::istringstream iss( "STATE_VECTOR min trpcage.list\nFITNESS vmin( trpcage ) + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << "test_daf_STATE_VECTOR_command_illegal_varname" << std::endl;
			//std::cout << e.msg() << std::endl;
			std::string expected_error_message = "Illegal name for variable, 'min' in the STATE_VECTOR command on line 1\nSTATE_VECTOR min trpcage.list";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}
	}

	void test_daf_STATE_VECTOR_command_function_varname() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );

		std::istringstream iss( "STATE_VECTOR sqrt trpcage.list\nFITNESS vmin( trpcage ) + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << "test_daf_STATE_VECTOR_command_function_varname" << std::endl;
			//std::cout << e.msg() << std::endl;
			std::string expected_error_message = //"Declaration of variable 'sqrt' in STATE_VECTOR conflicts with a function name.  Line 1\nSTATE_VECTOR sqrt trpcage.list";
				"Declaration of variable 'sqrt' in STATE_VECTOR command conflicts with a function name.  Line 1\n"
				"STATE_VECTOR sqrt trpcage.list";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}
	}

	void test_daf_STATE_VECTOR_command_duplicate_varname() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );

		std::istringstream iss( "STATE_VECTOR trpcage trpcage.list\nSTATE_VECTOR trpcage trpcage.list\nFITNESS vmin( trpcage ) + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << "test_daf_STATE_VECTOR_command_duplicate_varname" << std::endl;
			//std::cout << e.msg() << std::endl;
			std::string expected_error_message = //"Variable name trpcage appears multiple times in the DynamicAggregateFunction file.\nFirst occurrance was found on line 1.  Second occurrance found while reading STATE_VECTOR command\nSTATE_VECTOR trpcage trpcage.list\nLine # 2";
				"Variable name trpcage appears multiple times in the DynamicAggregateFunction file.\n"
				"Second occurrance found while reading a STATE_VECTOR command\n"
				"STATE_VECTOR trpcage trpcage.list\n"
				"Line # 2";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}

		}
	}

	void test_daf_POSE_ENERGY_command() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );
		daf->set_score_function( *sfxn );

		std::istringstream iss( "STATE_VECTOR trpcage trpcage.list\nPOSE_ENERGY trpcage_wt 1l2y.pdb\nFITNESS trpcage - trpcage_wt\n" );
		daf->initialize_from_input_file( ds, iss );

	}


	void test_daf_scalar_NPD_PROPERTY_command() {

		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );
		ds->add_npdpro_calculator_creator( new TestCalculatorCreator );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );

		std::istringstream iss( "STATE trpcage 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\nNPD_PROPERTY sasa_trpcage trpcage TESTCALC\nFITNESS trpcage + 2.0\n" );

		daf->initialize_from_input_file( ds, iss );

	}

	void test_daf_vector_NPD_PROPERTY_command() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );
		ds->add_npdpro_calculator_creator( new TestCalculatorCreator );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );


		std::istringstream iss( "STATE_VECTOR trpcage trpcage.list\nNPD_PROPERTY sasa_trpcage trpcage TESTCALC\nFITNESS vmin( trpcage ) + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( true );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cout << e.msg();
			TS_ASSERT( false );
		}
	}

	void test_daf_VECTOR_EXPRESSION_command() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );

		std::istringstream iss( "STATE_VECTOR trpcage trpcage.list\nVECTOR_EXPRESSION FOR x IN trpcage : trpcage_plus_2 = x + 2\nFITNESS vmin( trpcage_plus_2 ) + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( true );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cout << e.msg();
			TS_ASSERT( false );
		}
	}

	void test_daf_VECTOR_EXPRESSION_command_two_local_variables() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );

		std::istringstream iss( "STATE_VECTOR trpcage1 trpcage.list\nSTATE_VECTOR trpcage2 trpcage.list\nVECTOR_EXPRESSION FOR x IN trpcage1 , y IN trpcage2 : trpcage_plus_2 = x + y + 2\nFITNESS vmin( trpcage_plus_2 ) + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( true );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cout << e.msg();
			TS_ASSERT( false );
		}
	}

	void test_daf_VECTOR_EXPRESSION_command_two_local_variables_misplaced_comma() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );

		std::istringstream iss( "STATE_VECTOR trpcage1 trpcage.list\nSTATE_VECTOR trpcage2 trpcage.list\nVECTOR_EXPRESSION FOR x IN trpcage1, y IN trpcage2 : trpcage_plus_2 = x + y + 2\nFITNESS vmin( trpcage_plus_2 ) + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << "test_daf_VECTOR_EXPRESSION_command_two_local_variables_misplaced_comma" << std::endl;
			//std::cout << e.msg() << std::endl;
			std::string expected_error_message = "The vector-variable name 'trpcage1,'  given for local-variable 'x' does not belong to an already-declared vector variable.  Error in the VECTOR_EXPRESSION command on line 3\nVECTOR_EXPRESSION FOR x IN trpcage1, y IN trpcage2 : trpcage_plus_2 = x + y + 2";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}
	}


	void test_daf_VECTOR_EXPRESSION_command_nocap_IN() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );

		std::istringstream iss( "STATE_VECTOR trpcage trpcage.list\nVECTOR_EXPRESSION FOR x in trpcage : trpcage_plus_2 = x + 2\nFITNESS vmin( trpcage_plus_2 ) + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << e.msg();
			std::string expected_error_message = "Expected to read 'IN' in the DynamicAggregateFunction input file following the declaration of local variable 'x', but read 'in' in the VECTOR_EXPRESSION command on line 2\nVECTOR_EXPRESSION FOR x in trpcage : trpcage_plus_2 = x + 2";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}
	}

	void test_daf_VECTOR_EXPRESSION_command_duplicated_local_varname() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );

		std::istringstream iss( "STATE_VECTOR trpcage trpcage.list\nSTATE_VECTOR trpcage2 trpcage.list\nVECTOR_EXPRESSION FOR x IN trpcage , x IN trpcage2  : trpcage_plus_2 = x + 2\nFITNESS vmin( trpcage_plus_2 ) + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << "test_daf_VECTOR_EXPRESSION_command_duplicated_local_varname()" << std::endl;
			//std::cout << e.msg() << std::endl;
			std::string expected_error_message = "Local variable, 'x' in VECTOR_EXPRESSION command appears multiple times on line 3\nVECTOR_EXPRESSION FOR x IN trpcage , x IN trpcage2  : trpcage_plus_2 = x + 2";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}
	}

	void test_daf_VECTOR_VARIABLE_command() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );

		std::istringstream iss( "STATE trpcage1 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\nSTATE trpcage2 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\nVECTOR_VARIABLE trpcage_v = trpcage1 trpcage2\nVECTOR_EXPRESSION FOR x IN trpcage_v : trpcage_plus_2 = x + 2\nFITNESS vmin( trpcage_plus_2 ) + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( true );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cout << e.msg();
			TS_ASSERT( false );
		}
	}

	void test_daf_VECTOR_VARIABLE_command_bad_varname() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );

		std::istringstream iss( "STATE trpcage1 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\nSTATE trpcage2 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\nVECTOR_VARIABLE trpcage_v = trpcage1 trpcage3\nVECTOR_EXPRESSION FOR x IN trpcage_v : trpcage_plus_2 = x + 2\nFITNESS vmin( trpcage_plus_2 ) + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << "test_daf_VECTOR_VARIABLE_command_bad_varname" << std::endl;
			//std::cout << e.msg() << std::endl;
			std::string expected_error_message = "Unknown variable 'trpcage3' requested in the VECTOR_VARIABLE command on line 3\nVECTOR_VARIABLE trpcage_v = trpcage1 trpcage3";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}
	}

	void test_daf_VECTOR_VARIABLE_command_bad_vecvarname() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );

		std::istringstream iss( "STATE trpcage1 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\nSTATE trpcage2 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\nVECTOR_VARIABLE min = trpcage1 trpcage2\nVECTOR_EXPRESSION FOR x IN trpcage_v : trpcage_plus_2 = x + 2\nFITNESS vmin( trpcage_plus_2 ) + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << "test_daf_VECTOR_VARIABLE_command_bad_vecvarname" << std::endl;
			//std::cout << e.msg() << std::endl;
			std::string expected_error_message = "Illegal name for variable, 'min' in the VECTOR_VARIABLE command on line 3\nVECTOR_VARIABLE min = trpcage1 trpcage2";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}
	}

	void test_daf_VECTOR_VARIABLE_command_func_vecvarname() {
		DaemonSetOP ds = new DaemonSet;
		ScoreFunctionOP sfxn = getScoreFunction();
		ds->set_score_function( *sfxn );
		std::string corr_resfile_string( "2\nstart\n1 A PIKAA AP\n 2 A PIKAA FW EX ARO 1 LEVEL 4 EX ARO 2 LEVEL 4 EX_CUTOFF 1\n" );
		std::istringstream corr_resfile( corr_resfile_string );
		ds->set_entity_resfile( corr_resfile, "unnamed" );

		DynamicAggregateFunctionOP daf = new DynamicAggregateFunction;
		daf->add_file_contents( "1l2y.pdb", trp_cage_ideal() );
		daf->add_file_contents( "1l2y.correspondence.txt", "1 9 A\n2 13 A\n" );
		daf->add_file_contents( "1l2y.secondary.resfile", "NATRO\nstart\n9 A NATAA\n13 A NATAA\n" );
		daf->add_file_contents( "trpcage.list", "1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n" );

		std::istringstream iss( "STATE trpcage1 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n"
			"STATE trpcage2 1l2y.pdb 1l2y.correspondence.txt 1l2y.secondary.resfile\n"
			"VECTOR_VARIABLE exp = trpcage1 trpcage2\n"
			"VECTOR_EXPRESSION FOR x IN trpcage_v : trpcage_plus_2 = x + 2\n"
			"FITNESS vmin( trpcage_plus_2 ) + 2.0\n" );

		try {
			daf->initialize_from_input_file( ds, iss );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << "test_daf_VECTOR_VARIABLE_command_func_vecvarname" << std::endl;
			//std::cout << e.msg() << std::endl;
			std::string expected_error_message = //"Declaration of variable 'exp' in VECTOR_VARIABLE conflicts with a function name.  Line 3\nVECTOR_VARIABLE exp = trpcage1 trpcage2";
				"Declaration of variable 'exp' in VECTOR_VARIABLE command conflicts with a function name.  Line 3\n"
				"VECTOR_VARIABLE exp = trpcage1 trpcage2";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}
	}

	void test_EntityFunc_read_simple_file() {
		EntityFunc entfunc;
		entfunc.set_num_entity_elements( 1 );
		std::string efuncfile( "SCORE 1\n" );
		std::istringstream iss( efuncfile );
		entfunc.initialize_from_input_file( iss );

		Entity ent( "traits AA:1:F fitness 0.0" );

		core::Real score = entfunc.evaluate( ent );
		TS_ASSERT( score == 1.0 );

	}

	void test_EntityFunc_AA_SET_command() {
		EntityFunc entfunc;
		entfunc.set_num_entity_elements( 1 );
		std::string efuncfile( "AA_SET polar = { d,e,h,k,n,q,r,s,t}\nSCORE 1\n" );
		std::istringstream iss( efuncfile );
		entfunc.initialize_from_input_file( iss );

		Entity ent( "traits AA:1:F fitness 0.0" );

		core::Real score = entfunc.evaluate( ent );
		TS_ASSERT( score == 1.0 );

	}

	void test_EntityFunc_AA_SET_command_missing_commas() {
		EntityFunc entfunc;
		entfunc.set_num_entity_elements( 1 );
		std::string efuncfile( "AA_SET polar = { d e,h,k,n,q,r,s,t}\nSCORE 1\n" );
		std::istringstream iss( efuncfile );
		try {
			entfunc.initialize_from_input_file( iss );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << e.msg();
			std::string expected_error_message = "Expected to read comma  or a right curly brace, but found 'e'.\n"
				"Error encountered while reading AA_SET command\n"
				"AA_SET polar = { d e,h,k,n,q,r,s,t}\n"
				"Line # 1";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}

	}

	void test_EntityFunc_AA_SET_command_missing_equals_sign() {
		EntityFunc entfunc;
		entfunc.set_num_entity_elements( 1 );
		std::string efuncfile( "AA_SET polar { d, e,h,k,n,q,r,s,t}\nSCORE 1\n" );
		std::istringstream iss( efuncfile );
		try {
			entfunc.initialize_from_input_file( iss );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << e.msg();
			TS_ASSERT(
				"Expected to read an equals sign after reading amino-acid-set name'polar' but found '{'\n"
				"Error encountered while reading AA_SET command\n"
				"AA_SET polar { d, e,h,k,n,q,r,s,t}\n"
				"Line # 1"	== e.msg() );
		}

	}

	void test_EntityFunc_AA_SET_command_missing_right_curly_brace() {
		EntityFunc entfunc;
		entfunc.set_num_entity_elements( 1 );
		std::string efuncfile( "AA_SET polar = { d, e,h,k,n,q,r,s,t\nSCORE 1\n" );
		std::istringstream iss( efuncfile );
		try {
			entfunc.initialize_from_input_file( iss );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << e.msg();
			std::string expected_error_message =
				"Expected to read a right curly bracket ('}') or a 1-letter  amino acid code, but found an end-of-line.\n"
				"Error encountered while reading AA_SET command\n"
				"AA_SET polar = { d, e,h,k,n,q,r,s,t\n"
				"Line # 1";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}

	}

	void test_EntityFunc_SET_CONDITION_command_from_AA_SET() {
		EntityFunc entfunc;
		entfunc.set_num_entity_elements( 1 );
		std::string efuncfile( "AA_SET polar = { d,e,h,k,n,q,r,s,t}\nSET_CONDITION pol1 = ee_1 in polar\n SCORE pol1\n" );
		std::istringstream iss( efuncfile );
		entfunc.initialize_from_input_file( iss );

		Entity ent1( "traits AA:1:F fitness 0.0" );

		core::Real score = entfunc.evaluate( ent1 );
		TS_ASSERT( score == 0.0 );

		Entity ent2( "traits AA:1:R fitness 0.0" );

		score = entfunc.evaluate( ent2 );
		TS_ASSERT( score == 1.0 );

	}

	void test_EntityFunc_SET_CONDITION_command_from_aa_list() {
		EntityFunc entfunc;
		entfunc.set_num_entity_elements( 1 );
		std::string efuncfile( "SET_CONDITION pol1 = ee_1 in { d,e,h,k,n,q,r,s,t}\n SCORE pol1\n" );
		std::istringstream iss( efuncfile );
		entfunc.initialize_from_input_file( iss );

		Entity ent1( "traits AA:1:F fitness 0.0" );

		core::Real score = entfunc.evaluate( ent1 );
		TS_ASSERT( score == 0.0 );

		Entity ent2( "traits AA:1:R fitness 0.0" );

		score = entfunc.evaluate( ent2 );
		TS_ASSERT( score == 1.0 );

	}

	void test_EntityFunc_SET_CONDITION_command_from_AA_SET_bad_name() {
		EntityFunc entfunc;
		entfunc.set_num_entity_elements( 1 );
		std::string efuncfile( "AA_SET polar = { d,e,h,k,n,q,r,s,t}\nSET_CONDITION pol1 = ee_1 in ppolar\n SCORE pol1\n" );
		std::istringstream iss( efuncfile );
		try {
			entfunc.initialize_from_input_file( iss );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string expected_error_message = "Amino-acid-set name 'ppolar' has not previously been declared.\n"
				"Error encountered while reading SET_CONDITION command\n"
				"SET_CONDITION pol1 = ee_1 in ppolar\n"
				"Line # 2";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}


	}


	void test_EntityFunc_SET_CONDITION_command_from_aa_list_missing_comma() {
		EntityFunc entfunc;
		entfunc.set_num_entity_elements( 1 );
		std::string efuncfile( "SET_CONDITION pol1 = ee_1 in { d e,h,k,n,q,r,s,t}\n SCORE pol1\n" );
		std::istringstream iss( efuncfile );
		try {
			entfunc.initialize_from_input_file( iss );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << e.msg();
			std::string expected_error_message = "Expected to read comma  or a right curly brace, but found 'e'.\n"
				"Error encountered while reading SET_CONDITION command\n"
				"SET_CONDITION pol1 = ee_1 in { d e,h,k,n,q,r,s,t}\n"
				"Line # 1";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}
	}

	void test_EntityFunc_SET_CONDITION_command_from_aa_list_missing_lcurly() {
		EntityFunc entfunc;
		entfunc.set_num_entity_elements( 1 );
		std::string efuncfile( "SET_CONDITION pol1 = ee_1 in  d, e,h,k,n,q,r,s,t}\n SCORE pol1\n" );
		std::istringstream iss( efuncfile );
		try {
			entfunc.initialize_from_input_file( iss );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << e.msg();
			std::string expected_error_message = "Amino-acid-set name 'd,' has not previously been declared.\n"
				"Error encountered while reading SET_CONDITION command\n"
				"SET_CONDITION pol1 = ee_1 in  d, e,h,k,n,q,r,s,t}\n"
				"Line # 1";
			TS_ASSERT( expected_error_message == e.msg() );
			if ( expected_error_message != e.msg() ) {
				std::cout << "Actual error message\n\n" << e.msg() << std::endl;
			}
		}
	}

	void test_EntityFunc_SUB_EXPRESSION_command() {
		EntityFunc entfunc;
		entfunc.set_num_entity_elements( 1 );
		std::string efuncfile( "AA_SET polar = { d,e,h,k,n,q,r,s,t}\nSET_CONDITION pol1 = ee_1 in polar\nSUB_EXPRESSION scale_pol1 = -2 * pol1\nSCORE scale_pol1 - 4\n" );
		std::istringstream iss( efuncfile );
		entfunc.initialize_from_input_file( iss );

		Entity ent1( "traits AA:1:F fitness 0.0" );

		core::Real score = entfunc.evaluate( ent1 );
		TS_ASSERT( score == -4.0 );

		Entity ent2( "traits AA:1:R fitness 0.0" );

		score = entfunc.evaluate( ent2 );
		TS_ASSERT( score == -6.0 );
	}
};


