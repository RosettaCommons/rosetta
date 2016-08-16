// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/simulate_mpi.cxxtest.hh
/// @brief demonstate simulating mpi behavior in a single thread for unit testing
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <utility/mpi_util.hh>
#include <utility/SimulateMPI.hh>
#include <utility/excn/Exceptions.hh>

class SimulateMPITests : public CxxTest::TestSuite {

public:

	void setUp() {
		int nprocs = 4;
		utility::SimulateMPI::initialize_simulation(nprocs);
	}

	// @brief test default ctor
	void test_nprocs() {
		using namespace utility;
		TS_ASSERT(mpi_nprocs() == 4);
	}

	void test_string() {
		using namespace utility;
		// now we are rank 2
		SimulateMPI::set_mpi_rank(2);
		TS_ASSERT(mpi_rank() == 2);

		// send string to head node
		std::string orig_message("test_message");
		send_string_to_node(0, orig_message);

		SimulateMPI::set_mpi_rank(0);
		std::string message(receive_string_from_node(2));
		TS_ASSERT(message == orig_message);
		try {
			std::string m2 = receive_string_from_node(2);
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			TS_ASSERT_EQUALS( e.msg(), "Could not retrieve a smpi_string message to 0 from 2; message queue is empty" );
		}
	}

	void test_char() {
		using namespace utility;
		SimulateMPI::set_mpi_rank(2);
		TS_ASSERT(mpi_rank() == 2);

		char orig_message('a');
		send_char_to_node(0, orig_message);

		SimulateMPI::set_mpi_rank(0);
		char message(receive_char_from_node(2));
		TS_ASSERT(message == orig_message);
		try {
			/*char m2 =*/ receive_char_from_node(2);
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			TS_ASSERT_EQUALS( e.msg(), "Could not retrieve a smpi_char message to 0 from 2; message queue is empty" );
		}
	}

	void test_integer() {
		using namespace utility;
		SimulateMPI::set_mpi_rank(2);
		TS_ASSERT(mpi_rank() == 2);

		int orig_message(10);
		send_integer_to_node(0, orig_message);

		SimulateMPI::set_mpi_rank(0);
		int message(receive_integer_from_node(2));
		TS_ASSERT(message == orig_message);

		try {
			/*int m2 =*/ receive_integer_from_node(2);
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			TS_ASSERT_EQUALS( e.msg(), "Could not retrieve a smpi_integer message to 0 from 2; message queue is empty" );
		}
	}

	void test_integers() {
		using namespace utility;
		SimulateMPI::set_mpi_rank(2);
		TS_ASSERT(mpi_rank() == 2);

		vector1< int > orig_message;
		orig_message.push_back(10);
		send_integers_to_node(0, orig_message);

		SimulateMPI::set_mpi_rank(0);
		vector1< int > message(receive_integers_from_node(2));
		TS_ASSERT(message[1] == orig_message[1]);

		try {
			vector1< int > m2 = receive_integers_from_node(2);
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			TS_ASSERT_EQUALS( e.msg(), "Could not retrieve a smpi_integers message to 0 from 2; message queue is empty" );
		}
	}


	void test_double() {
		using namespace utility;
		SimulateMPI::set_mpi_rank(2);
		TS_ASSERT(mpi_rank() == 2);

		double orig_message(10);
		send_double_to_node(0, orig_message);

		SimulateMPI::set_mpi_rank(0);
		double message(receive_double_from_node(2));
		TS_ASSERT(message == orig_message);
		try {
			/*double m2 =*/ receive_double_from_node(2);
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			TS_ASSERT_EQUALS( e.msg(), "Could not retrieve a smpi_double message to 0 from 2; message queue is empty" );
		}

	}

	void test_doubles() {
		using namespace utility;
		SimulateMPI::set_mpi_rank(2);
		TS_ASSERT(mpi_rank() == 2);

		vector1< double > orig_message;
		orig_message.push_back(10);
		send_doubles_to_node(0, orig_message);

		SimulateMPI::set_mpi_rank(0);
		vector1< double > message(receive_doubles_from_node(2));
		TS_ASSERT(message[1] == orig_message[1]);
		try {
			vector1< double > m2 = receive_doubles_from_node(2);
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			TS_ASSERT_EQUALS( e.msg(), "Could not retrieve a smpi_doubles message to 0 from 2; message queue is empty" );
		}
	}


};
