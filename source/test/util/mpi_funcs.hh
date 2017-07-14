// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/util/mpi_funcs.hh
/// @brief  common utility functions for testing MPI functionality
/// @author Andrew Leaver-Fay
/// @author Andy Watkins

#ifndef INCLUDED_util_mpi_funcs_HH
#define INCLUDED_util_mpi_funcs_HH

#include <string>
#include <utility/SimulateMPI.hh>
#include <utility/mpi_util.hh>
#include <utility/excn/Exceptions.hh>

void ts_assert_mpi_buffer_has_string(
	core::Size source,
	std::string message_tag,
	std::string expected_message
);


std::string ts_assert_mpi_buffer_has_string(
	core::Size source,
	std::string message_tag
);

void ts_assert_mpi_buffer_has_integer(
	core::Size source,
	std::string message_tag,
	int expected_message
);

void ts_assert_mpi_buffer_has_size(
	core::Size source,
	std::string message_tag,
	core::Size expected_message
);

double ts_assert_mpi_buffer_has_double(
	core::Size source,
	std::string message_tag
);

void ts_assert_mpi_buffer_has_string(
	core::Size source,
	std::string message_tag,
	std::string expected_message
) {
	try {
		std::string msg = utility::receive_string_from_node( source );
		TS_ASSERT( msg == expected_message );
		if ( msg != expected_message ) {
			std::cerr << "SimulateMPI string for tag \"" << message_tag << "\" did not match expected string:\n";
			std::cerr << "Expected: " << expected_message << "\n";
			std::cerr << "Actual: " << msg << "\n";
		}
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::cerr << "Exception caught for tag \"" << message_tag << "\": " << e.msg() << std::endl;
		TS_ASSERT( false );
	}
}

std::string ts_assert_mpi_buffer_has_string(
	core::Size source,
	std::string message_tag
) {
	std::string msg;
	try {
		msg = utility::receive_string_from_node( source );
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::cerr << "Exception caught for tag \"" << message_tag << "\": " << e.msg() << std::endl;
		TS_ASSERT( false );
	}
	return msg;
}

void ts_assert_mpi_buffer_has_integer(
	core::Size source,
	std::string message_tag,
	int expected_message
) {
	try {

		//std::cout << "Testing if we receive " << message_tag << expected_message << " from " << source << std::endl;
		int msg = utility::receive_integer_from_node( source );
		//std::cout << "message was " << msg << std::endl;
		TS_ASSERT( msg == expected_message );
		if ( msg != expected_message ) {
			std::cerr << "SimulateMPI int for tag \"" << message_tag << "\" did not match expected int:\n";
			std::cerr << "Expected: \"" << expected_message << "\"\n";
			std::cerr << "Actual: \"" << msg << "\"\n";
		}
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::cerr << "Exception caught for tag \"" << message_tag << "\": " << e.msg() << std::endl;
		TS_ASSERT( false );
	}
}

void ts_assert_mpi_buffer_has_size(
	core::Size source,
	std::string message_tag,
	core::Size expected_message
) {
	try {

		//std::cout << "Testing if we receive " << message_tag << expected_message << " from " << source << std::endl;
		core::Size msg = utility::receive_size_from_node( source );
		//std::cout << "message was " << msg << std::endl;
		TS_ASSERT( msg == expected_message );
		if ( msg != expected_message ) {
			std::cerr << "SimulateMPI string for tag \"" << message_tag << "\" did not match expected string:\n";
			std::cerr << "Expected: " << expected_message << "\n";
			std::cerr << "Actual: " << msg << "\n";
		}
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::cerr << "Exception caught for tag \"" << message_tag << "\": " << e.msg() << std::endl;
		TS_ASSERT( false );
	}
}

double ts_assert_mpi_buffer_has_double(
	core::Size source,
	std::string message_tag
) {
	double msg( 0 );
	try {
		msg = utility::receive_double_from_node( source );
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::cerr << "Exception caught for tag \"" << message_tag << "\": " << e.msg() << std::endl;
		TS_ASSERT( false );
	}
	return msg;
}

void ts_assert_mpi_buffer_has_integers(
	core::Size source,
	std::string message_tag,
	utility::vector1< int > const & expected_message
) {
	try {

		//std::cout << "Testing if we receive " << message_tag << expected_message << " from " << source << std::endl;
		utility::vector1< int > msg = utility::receive_integers_from_node( source );
		//std::cout << "message was " << msg << std::endl;
		TS_ASSERT( msg == expected_message );
		if ( msg != expected_message ) {
			std::cerr << "SimulateMPI integer vector for tag \"" << message_tag << "\" did not match expected integer vector:\n";
			std::cerr << "Expected: [ ";
			for ( core::Size ii = 1; ii <= expected_message.size(); ++ii ) {
				if ( ii != 1 ) std::cerr << ", ";
				std::cerr << expected_message[ ii ];
			}
			std::cerr << " ]\n";
			std::cerr << "Actual: [ ";
			for ( core::Size ii = 1; ii <= msg.size(); ++ii ) {
				if ( ii != 1 ) std::cerr << ", ";
				std::cerr << msg[ ii ];
			}
			std::cerr << " ]\n";
		}
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::cerr << "Exception caught for tag \"" << message_tag << "\": " << e.msg() << std::endl;
		TS_ASSERT( false );
	}
}

void ts_assert_mpi_buffer_has_sizes(
	core::Size source,
	std::string message_tag,
	utility::vector1< core::Size > const & expected_message
) {
	try {

		//std::cout << "Testing if we receive " << message_tag << expected_message << " from " << source << std::endl;
		utility::vector1< core::Size > msg = utility::receive_sizes_from_node( source );
		//std::cout << "message was " << msg << std::endl;
		TS_ASSERT( msg == expected_message );
		if ( msg != expected_message ) {
			std::cerr << "SimulateMPI size vector for tag \"" << message_tag << "\" did not match expected size vector:\n";
			std::cerr << "Expected: [ ";
			for ( core::Size ii = 1; ii <= expected_message.size(); ++ii ) {
				if ( ii != 1 ) std::cerr << ", ";
				std::cerr << expected_message[ ii ];
			}
			std::cerr << " ]\n";
			std::cerr << "Actual: [ ";
			for ( core::Size ii = 1; ii <= msg.size(); ++ii ) {
				if ( ii != 1 ) std::cerr << ", ";
				std::cerr << msg[ ii ];
			}
			std::cerr << " ]\n";
		}
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::cerr << "Exception caught for tag \"" << message_tag << "\": " << e.msg() << std::endl;
		TS_ASSERT( false );
	}
}

#endif
