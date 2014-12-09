// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/devel/mmt_msd/MMTDriver.cxxtest.hh
/// @brief  test suite for devel::mmt_msd::MMTDriver using the SimultateMPI utilities
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Unit headers
#include <devel/mmt_msd/MMTDriver.hh>
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

// Protocols headers
#include <protocols/genetic_algorithm/Entity.hh>
#include <protocols/genetic_algorithm/EntityRandomizer.hh>

// Utility headers
#include <utility/SimulateMPI.hh>
#include <utility/mpi_util.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/io/FileContentsMap.hh>

// --------------- Test Class --------------- //

// return entities in a perscribed order
class BogusRandomizer : public protocols::genetic_algorithm::PositionSpecificRandomizer
{
public:
	typedef protocols::genetic_algorithm::Entity Entity;
	typedef protocols::genetic_algorithm::EntityOP EntityOP;

public:

	virtual ~BogusRandomizer() {}

	std::string
	eestring( std::string seq ) {
		std::ostringstream ostr;
		ostr << "traits";
		for ( core::Size ii = 0; ii < seq.length(); ++ii ) {
			ostr << " AA:" << ii+1 << ":" << seq[ii];
		}
		ostr << " fitness 0.0";
		return ostr.str();
	};

	EntityOP next_entity() {
		std::string seq = seqs_.front();
		seqs_.pop_front();
		EntityOP temp( new Entity(eestring( seq )) );
		EntityOP retval( new Entity );
		retval->set_traits( temp->traits() );
		return retval;
	}

	virtual EntityOP random_entity() {
		return next_entity();
	}

	virtual void mutate( protocols::genetic_algorithm::Entity & entity )  {
		EntityOP next = next_entity();
		entity.set_traits( next->traits() );
	}

	virtual core::Size library_size() const { return seqs_.size(); }

	void append_sequence( std::string const & newseq ) {
		seqs_.push_back( newseq );
	}



private:
	std::list< std::string > seqs_;

};

typedef utility::pointer::shared_ptr< BogusRandomizer > BogusRandomizerOP;

class MMTDriverTests : public CxxTest::TestSuite {
public:

	typedef core::Size   Size;
	typedef core::Real   Real;
	typedef std::list< std::pair< std::string, core::Real > > SeqsAndEnergies;

private:
	utility::io::FileContentsMapOP fc_;

public:

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

	void setup_driver( devel::mmt_msd::MMTDriver & mmt_driver )
	{
		mmt_driver.set_nworkers( 1 );
		mmt_driver.set_ngenerations( 3 );
		mmt_driver.set_pop_size( 2 );
		mmt_driver.set_frac_by_recomb( 0 );

		BogusRandomizerOP rand( new BogusRandomizer );
		rand->append_sequence( "AG" ); // gen1a
		rand->append_sequence( "AQ" ); // gen1b
		rand->append_sequence( "MQ" ); // gen2
		rand->append_sequence( "MG" ); // gen3
		mmt_driver.set_randomizer( rand );

		fc_ = utility::io::FileContentsMapOP( new utility::io::FileContentsMap );
		fc_->set_file_contents( "ubq.pdb", ubq_twores_string() );
		fc_->set_file_contents( "ubq.corr", "1 1 _\n2 2 _\n" );
		fc_->set_file_contents( "ubq_entity.resfile", "2\nstart\n1 A PIKAA AM\n2 A PIKAA GQ\n" );
		fc_->set_file_contents( "ubq.2res", "NATRO\nstart\n" );
		fc_->set_file_contents( "ubq.daf", "STATE ubq ubq.pdb ubq.corr ubq.2res\nFITNESS ubq\n");
		mmt_driver.set_file_contents( fc_ );
		mmt_driver.set_n_results_to_output( 1 );
		mmt_driver.set_daf_fname( "ubq.daf" );
		mmt_driver.set_entity_resfile_fname( "ubq_entity.resfile" );
	}

	void initialize_handshake_mpi_messages() {
		utility::SimulateMPI::set_mpi_rank( 1 );
		utility::send_integer_to_node( 0, devel::mmt_msd::handshake_acknowledged );
		utility::send_integer_to_node( 0, 2 );
	}

	void queue_ready_for_work_mpi_messages() {
		utility::SimulateMPI::set_mpi_rank( 1 );
		utility::send_integer_to_node( 0, 1 );
		utility::send_integer_to_node( 0, devel::mmt_msd::requesting_new_job );
	}

	void queue_job_completed_mpi_messages( std::string seq, core::Real energy ) {
		utility::SimulateMPI::set_mpi_rank( 1 );
		utility::send_integer_to_node( 0, 1 );
		utility::send_integer_to_node( 0, devel::mmt_msd::job_complete );
		utility::send_integer_to_node( 0, 1 ); // state index
		utility::send_string_to_node( 0, seq );
		utility::send_double_to_node( 0, energy );
		utility::send_double_to_node( 0, 0.125 ); // running time
		utility::send_integer_to_node( 0, 0 ); // num npd properties
	}

	//void queue_repack_ubqtwores_mpi_messages() {
	//	utility::SimulateMPI::set_mpi_rank( 0 );
	//
	//	// corresponding to the series of communications in MMTDriver::receive_state_input_data
	//	utility::send_integer_to_node( 1, 1 ); // state index
	//	utility::send_string_to_node( 1, "1ubq.pdb" );
	//	utility::send_string_to_node( 1, ubq_twores_string() );
	//	utility::send_string_to_node( 1, "bogus.corr" );
	//	utility::send_string_to_node( 1, "1 1 _\n" );
	//	utility::send_string_to_node( 1, "bogus.2res" );
	//	utility::send_string_to_node( 1, "NATRO\nstart\n" );
	//	utility::send_integer_to_node( 1, 0 ); // no npd properties
	//}


	//void queue_request_to_save_lastgen_sequence( std::string const & seq ) {
	//	utility::send_integer_to_node( 1, devel::mmt_msd::result_from_last_generation_needs_saving );
	//	utility::send_integer_to_node( 1, 1 );
	//	utility::send_string_to_node( 1, seq );
	//}

	void initialize_main_opt_loop_messages(
		SeqsAndEnergies const & gen_seqs_and_energies
	) {
		queue_ready_for_work_mpi_messages();
		for ( SeqsAndEnergies::const_iterator
				iter = gen_seqs_and_energies.begin(),
				iter_end = gen_seqs_and_energies.end();
				iter != iter_end; ++iter ) {
			queue_ready_for_work_mpi_messages();
		}

		for ( SeqsAndEnergies::const_iterator
				iter = gen_seqs_and_energies.begin(),
				iter_end = gen_seqs_and_energies.end();
				iter != iter_end; ++iter ) {
			queue_job_completed_mpi_messages( iter->first, iter->second );
		}
	}

	void initialize_retrieve_optimal_pdb_messages( std::string const & )
	{
		utility::send_integer_to_node( 0, devel::mmt_msd::recovery_successful );
		utility::send_string_to_node( 0, fc_->get_file_contents_ref( "ubq.pdb" ) );
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
				std::cerr << "SimulateMPI string for tag '" << message_tag << "' did not match expected string:";
				std::cerr << "Expected: " << expected_message << "\n";
				std::cerr << "Actual: " << msg << "\n";
			}
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cerr << "Exception caught for tag '" << message_tag << "': " << e.msg() << std::endl;
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
			std::cerr << "Exception caught for tag '" << message_tag << "': " << e.msg() << std::endl;
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
				std::cerr << "SimulateMPI integer for tag '" << message_tag << "' did not match expected integer:";
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
		double msg;
		try {
			msg = utility::receive_double_from_node( source );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cerr << "Exception caught for tag " << message_tag << ": " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		return msg;
	}

	void verify_ubq_state_info( std::string const & event_prefix )
	{
		ts_assert_mpi_buffer_has_integer( 0, event_prefix + "state index message", 1 );
		ts_assert_mpi_buffer_has_string( 0, event_prefix + "ubq pdb name", "ubq.pdb" );
		ts_assert_mpi_buffer_has_string( 0, event_prefix + "ubq pdb string", fc_->get_file_contents_ref( "ubq.pdb" ) );
		ts_assert_mpi_buffer_has_string( 0, event_prefix + "ubq corr fname", "ubq.corr" );
		ts_assert_mpi_buffer_has_string( 0, event_prefix + "ubq corr string", fc_->get_file_contents_ref( "ubq.corr" ) );
		ts_assert_mpi_buffer_has_string( 0, event_prefix + "ubq corr fname", "ubq.2res" );
		ts_assert_mpi_buffer_has_string( 0, event_prefix + "ubq corr string", fc_->get_file_contents_ref( "ubq.2res" ) );
		ts_assert_mpi_buffer_has_integer( 0, event_prefix + "num npd properties", 0 );
	}

	void verify_one_generation_mpi_messages(
		core::Size gen_index,
		SeqsAndEnergies const & seqs_and_energies
	)
	{
		utility::SimulateMPI::set_mpi_rank( 1 );
		std::string genname = "generation " + utility::to_string( gen_index ) + " ";
		ts_assert_mpi_buffer_has_integer( 0, genname + "new generation message", devel::mmt_msd::new_generation );
		for ( SeqsAndEnergies::const_iterator
				iter = seqs_and_energies.begin(),
				iter_end = seqs_and_energies.end();
				iter != iter_end; ++iter ) {
			ts_assert_mpi_buffer_has_integer( 0, genname + "job ready message", devel::mmt_msd::new_job_ready );
			verify_ubq_state_info( genname );
			ts_assert_mpi_buffer_has_string( 0, genname + "sequence " + iter->first, iter->first );
		}
		ts_assert_mpi_buffer_has_integer( 0, genname + "all jobs assigned", devel::mmt_msd::generation_complete );
	}

	void verify_keep_sequence( std::string const & seq ) {
		ts_assert_mpi_buffer_has_integer( 0, "keep sequence message for " + seq + " part 1", devel::mmt_msd::result_from_last_generation_needs_saving );
		ts_assert_mpi_buffer_has_integer( 0, "keep sequence message for " + seq + " part 2", 1 );
		ts_assert_mpi_buffer_has_string( 0, "keep sequence message for " + seq + " part 3", seq );
	}

	void verify_discard_old_sequence( std::string const & seq ) {
		ts_assert_mpi_buffer_has_integer( 0, "discard sequence message for " + seq + " part 1", devel::mmt_msd::old_result_can_be_discarded );
		ts_assert_mpi_buffer_has_integer( 0, "discard sequence message for " + seq + " part 2", 1 );
		ts_assert_mpi_buffer_has_string( 0, "discard sequence message for " + seq + " part 3", seq );
	}

	void test_MMTDriver_end_to_end() {
		using namespace devel::mmt_msd;

		MMTDriver mmt_driver;
		setup_driver( mmt_driver );

		initialize_handshake_mpi_messages();

		SeqsAndEnergies gen1;
		gen1.push_back( std::make_pair( "AG", 1.0 ));
		gen1.push_back( std::make_pair( "AQ", 2.0 ));
		initialize_main_opt_loop_messages( gen1 );

		SeqsAndEnergies gen2;
		gen2.push_back( std::make_pair( "MQ", -1.0 ));
		initialize_main_opt_loop_messages( gen2 );

		SeqsAndEnergies gen3;
		gen3.push_back( std::make_pair( "MG", 0.0 ));
		initialize_main_opt_loop_messages( gen3 );

		initialize_retrieve_optimal_pdb_messages( "MQ" );

		utility::SimulateMPI::set_mpi_rank( 0 );
		std::map< std::string, std::string > solutions;
		try {
			mmt_driver.setup();
			mmt_driver.run();
			solutions = mmt_driver.retrieve_optimal_solutions();

		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cerr << e.msg() << std::endl;
			TS_ASSERT( false );
		}

		/// OK! let's see whether it worked.

		/// 1st handshake
		utility::SimulateMPI::set_mpi_rank( 1 );
		ts_assert_mpi_buffer_has_integer( 0, "test_mmtdriver_end_to_end handshake 1", devel::mmt_msd::handshake_begin );
		ts_assert_mpi_buffer_has_string(  0, "test_mmtdriver_end_to_end handshake 2", "ubq_entity.resfile" );
		ts_assert_mpi_buffer_has_string(  0, "test_mmtdriver_end_to_end handshake 3", fc_->get_file_contents_ref( "ubq_entity.resfile" ) );

		verify_one_generation_mpi_messages( 1, gen1 );
		verify_keep_sequence( "AG" );

		verify_one_generation_mpi_messages( 2, gen2 );
		verify_keep_sequence( "MQ" );
		verify_discard_old_sequence( "AG" );

		verify_one_generation_mpi_messages( 3, gen3 );

		// ok, now we should be asking for the optimal sequence, MQ
		ts_assert_mpi_buffer_has_integer( 0, "recover MQ pose initial request", devel::mmt_msd::recover_pose_for_result );
		verify_ubq_state_info( "recover MQ pose" );
		ts_assert_mpi_buffer_has_string( 0, "recover pose, sequence", "MQ" );

		// and finally, verify the spin-down mesage
		ts_assert_mpi_buffer_has_integer( 0, "spin down", devel::mmt_msd::spin_down );


	}


};
