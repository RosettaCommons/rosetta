// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/devel/mmt_msd/MMTReceiver.cxxtest.hh
/// @brief  test suite for devel::mmt_msd::MMTReceiver using the SimultateMPI utilities
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Unit headers
#include <devel/mmt_msd/MMTReceiver.hh>

/// Project headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/import_pose_options.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/scmin/SidechainStateAssignment.hh>

// Utility headers
#include <utility/SimulateMPI.hh>
#include <utility/mpi_util.hh>
#include <utility/excn/Exceptions.hh>

// --------------- Test Class --------------- //

class MMTReceiverTests : public CxxTest::TestSuite {

public:

	typedef core::Size   Size;
	typedef core::Real   Real;


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		core_init();
		utility::SimulateMPI::initialize_simulation( 2 );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void initialize_handshake_mpi_messages() {
		utility::SimulateMPI::set_mpi_rank( 0 );
		utility::send_integer_to_node( 1, devel::mmt_msd::handshake_begin );
		utility::send_string_to_node( 1, "bogus.resfile" );
		utility::send_string_to_node( 1, "1\nPIKAA AG\nstart\n" );
	}

	void queue_repack_ubqtwores_mpi_messages() {
		utility::SimulateMPI::set_mpi_rank( 0 );

		// corresponding to the series of communications in MMTReceiver::receive_state_input_data
		utility::send_integer_to_node( 1, 1 ); // state index
		utility::send_string_to_node( 1, "1ubq.pdb" );
		utility::send_string_to_node( 1, ubq_twores_string() );
		utility::send_string_to_node( 1, "bogus.corr" );
		utility::send_string_to_node( 1, "1 1 _\n" );
		utility::send_string_to_node( 1, "bogus.2res" );
		utility::send_string_to_node( 1, "NATRO\nstart\n" );
		utility::send_integer_to_node( 1, 0 ); // no npd properties
	}

	void queue_one_sequence_generation( std::string const & seq ) {
		utility::SimulateMPI::set_mpi_rank( 0 );
		utility::send_integer_to_node( 1, devel::mmt_msd::new_generation );
		utility::send_integer_to_node( 1, devel::mmt_msd::new_job_ready );

		queue_repack_ubqtwores_mpi_messages();
		// send the sequence
		utility::send_string_to_node( 1, seq );

		// ok, the job will end, and the MMTReceiver will want to send back results
		// and the MMTReceiver needs no signal from us to go through this process

		// now send the generation complete message
		utility::send_integer_to_node( 1, devel::mmt_msd::generation_complete );

	}

	void queue_request_to_save_lastgen_sequence( std::string const & seq ) {
		utility::send_integer_to_node( 1, devel::mmt_msd::result_from_last_generation_needs_saving );
		utility::send_integer_to_node( 1, 1 );
		utility::send_string_to_node( 1, seq );
	}

	void initialize_main_opt_loop_messages() {

		queue_one_sequence_generation( "G" );
		// now request that the node save that sequence
		queue_request_to_save_lastgen_sequence( "G" );

		queue_one_sequence_generation( "A" );
		// now request that the node save sequence "A"
		queue_request_to_save_lastgen_sequence( "A" );
		// and then discard sequence "G"
		utility::send_integer_to_node( 1, devel::mmt_msd::old_result_can_be_discarded );
		utility::send_integer_to_node( 1, 1 );
		utility::send_string_to_node( 1, "G" );

		// now retreive the pdb for the "A" sequence
		utility::send_integer_to_node( 1, devel::mmt_msd::recover_pose_for_result );
		queue_repack_ubqtwores_mpi_messages();
		utility::send_string_to_node( 1, "A" );

		// and finally send the spin-down message
		utility::send_integer_to_node( 1, devel::mmt_msd::spin_down );
	}


	void ts_assert_mpi_buffer_has_string(
		core::Size source,
		std::string message_tag,
		std::string expected_message
	)
	{
		try {
			std::string msg = utility::receive_string_from_node( source );
			TS_ASSERT( msg == expected_message );
			if ( msg != expected_message ) {
				std::cerr << "SimulateMPI string for tag " << message_tag << " did not match expected string:";
				std::cerr << "Expected: " << expected_message << "\n";
				std::cerr << "Actual: " << msg << "\n";
			}
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cerr << "Exception caught for tag " << message_tag << ": " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
	}

	std::string ts_assert_mpi_buffer_has_string(
		core::Size source,
		std::string message_tag
	)
	{
		std::string msg;
		try {
			msg = utility::receive_string_from_node( source );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cerr << "Exception caught for tag " << message_tag << ": " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		return msg;
	}

	void ts_assert_mpi_buffer_has_integer(
		core::Size source,
		std::string message_tag,
		int expected_message
	)
	{
		try {
			int msg = utility::receive_integer_from_node( source );
			TS_ASSERT( msg == expected_message );
			if ( msg != expected_message ) {
				std::cerr << "SimulateMPI string for tag " << message_tag << " did not match expected string:";
				std::cerr << "Expected: " << expected_message << "\n";
				std::cerr << "Actual: " << msg << "\n";
			}
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cerr << "Exception caught for tag " << message_tag << ": " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
	}

	double ts_assert_mpi_buffer_has_double(
		core::Size source,
		std::string message_tag
	)
	{
		double msg( 0 );
		try {
			msg = utility::receive_double_from_node( source );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cerr << "Exception caught for tag " << message_tag << ": " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		return msg;
	}

	void test_provide_an_empty_test() {
		TS_ASSERT( true );
	}

	void dont_test_MMTReceiver_handshake() {
		using namespace devel::mmt_msd;

		initialize_handshake_mpi_messages();
		utility::SimulateMPI::set_mpi_rank( 1 );
		MMTReceiver mmtr;
		mmtr.set_max_capacity( 1 );
		try {
			mmtr.initial_handshake();
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cerr << e.msg() << std::endl;
			TS_ASSERT( false );
		}

// ok, so the mmtr should've sent node 0 a few messages
		utility::SimulateMPI::set_mpi_rank( 0 );
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_handshake::handshake acknowledgement", devel::mmt_msd::handshake_acknowledged );
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_handshake::handshake max capacity", 1 );
		std::cout << "test_MMTReceiver_handshake complete" << std::endl;
	}

	void dont_test_MMTReceiver_end_to_end() {
		using namespace devel::mmt_msd;

		initialize_handshake_mpi_messages();
		initialize_main_opt_loop_messages();

		utility::SimulateMPI::set_mpi_rank( 1 );

		utility::SimulateMPI::set_mpi_rank( 1 );
		MMTReceiver mmtr;
		mmtr.set_max_capacity( 1 );
		try {
			mmtr.initial_handshake();
			mmtr.main_optimization_loop();

		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cerr << e.msg() << std::endl;
			TS_ASSERT( false );
		}

/// OK! let's see whether it worked.

/// handshake messages:
		utility::SimulateMPI::set_mpi_rank( 0 );
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::handshake acknowledgement", devel::mmt_msd::handshake_acknowledged );
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::handshake max capacity", 1 );

		/// ok, generation 1
		/// 1. node 1 should've identified itself as ready for a new job by sending its index
		/// 2. node 1 should've finished, and shipped back energies and npd properties
		/// 3. node 1 again should have identified itself and we say "generation ended"
		/// then we say "keep seq 'G'" and node 1 says nothing
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::gen1 client contact #1", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::gen1 work request #1", devel::mmt_msd::requesting_new_job );

		// energy, running time, and npd properties for sequence "G"
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::gen1 client contact #2", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::gen1 job complete #1", devel::mmt_msd::job_complete );
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::gen1 state index", 1);
		ts_assert_mpi_buffer_has_string( 1, "test_MMTReceiver_end_to_end::gen1 sequence", "G");
		ts_assert_mpi_buffer_has_double( 1, "test_MMTReceiver_end_to_end::gen1 energy" );
		double gen1_seq1_running_time = ts_assert_mpi_buffer_has_double( 1, "test_MMTReceiver_end_to_end::gen1 running time" );
		TS_ASSERT( gen1_seq1_running_time > 0 );
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::gen1 num_npd_props", 0);

		// ok, done sending packing results; node 1 should request a new job here.
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::gen1 client contact #3", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::gen1 work request #2", devel::mmt_msd::requesting_new_job );


		// ok, generation 2, again, same progression as in generation 1, but
		// expecting the results to correspond to the "A" sequence
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::gen2 client contact #1", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::gen2 work request #1", devel::mmt_msd::requesting_new_job );

		// energy, running time, and npd properties
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::gen2 client contact #2", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::gen2 job complete #1", devel::mmt_msd::job_complete );
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::gen2 state index", 1);
		ts_assert_mpi_buffer_has_string( 1, "test_MMTReceiver_end_to_end::gen2 sequence", "A");
		ts_assert_mpi_buffer_has_double( 1, "test_MMTReceiver_end_to_end::gen2 energy" );
		double gen2_seq1_running_time = ts_assert_mpi_buffer_has_double( 1, "test_MMTReceiver_end_to_end::gen2 running time" );
		TS_ASSERT( gen2_seq1_running_time > 0 );
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::gen2 num_npd_props", 0);


		// ok, done sending packing results; node 1 should request a new job here.
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::gen2 client contact #3", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end::gen2 work request #2", devel::mmt_msd::requesting_new_job );

		// ok, generation 2 ends, and now we say keep 'A' and discard 'G'
		// again, node 1 says nothing

		// up next: node 1 will ship us a pdb per our request
		ts_assert_mpi_buffer_has_integer( 1, "test_MMTReceiver_end_to_end pose recovery for sequence 'A'", devel::mmt_msd::recovery_successful );
		std::string pdb_string = ts_assert_mpi_buffer_has_string( 1, "test_MMTReceiver_end_to_end sequence 'A' pdb" );

		core::pose::Pose A_pose;
		core::import_pose::ImportPoseOptions import_opts;
		core::import_pose::pose_from_pdbstring( A_pose, pdb_string, import_opts, "test_MMTReceiver_end_to_end_A.pdb" );

		TS_ASSERT_EQUALS( A_pose.total_residue(), 2 );
		TS_ASSERT_EQUALS( A_pose.residue(1).aa(), core::chemical::aa_ala );
		TS_ASSERT_EQUALS( A_pose.residue(2).aa(), core::chemical::aa_gln );

		std::cout << "test_MMTReceiver_end_to_end complete" << std::endl;

	}


};
