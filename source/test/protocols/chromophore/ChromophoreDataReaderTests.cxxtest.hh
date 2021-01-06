// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/chromophore/ChromophoreDataReaderTests.cxxtest.hh
/// @brief  unit tests for ChromophoreDataReader class
/// @author Nina Bozhanova (nbozhanova@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/chromophore/ChromophoreDataReader.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("ChromophoreDataReaderTests");


class ChromophoreDataReaderTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		cro_data_reader_ = protocols::chromophore::ChromophoreDataReaderOP (new protocols::chromophore::ChromophoreDataReader());
	}

	void tearDown(){
		cro_data_reader_ = 0;
	}

	void test_nothing_to_parse_empty() {
		std::string empty_line("");
		std::istringstream test_stream(empty_line);
		try {
			// Mute BACKTRACE message
			set_throw_on_next_assertion_failure();
			// This line should throw the exception if everything works correctly
			cro_data_reader_->create_residue_names_maps(test_stream);
			// This line will fail the test if the exception is not thrown
			TS_ASSERT (false);
		}
catch (utility::excn::Exception & e ) {
	// Expected error message
	std::string const error = "No information in the input file. Please check your input file.";
	TS_ASSERT_STRING_CONTAINS(e.msg(), error);
}
	}

	void test_nothing_to_parse_comment_line() {
		std::string comment_line("#N1   N  ");
		std::istringstream test_stream(comment_line);
		try {
			// Mute BACKTRACE message
			set_throw_on_next_assertion_failure();
			// This line should throw the exception if everything works correctly
			cro_data_reader_->create_residue_names_maps(test_stream);
			// This line will fail the test if the exception is not thrown
			TS_ASSERT (false);
		}
catch (utility::excn::Exception & e ) {
	// Expected error message
	std::string const error = "No information in the input file. Please check your input file.";
	TS_ASSERT_STRING_CONTAINS(e.msg(), error);
}

	}

	void test_parse_correct_line() {
		std::string meaningful_line("N1   N    1");
		std::istringstream test_stream(meaningful_line);
		cro_data_reader_->create_residue_names_maps(test_stream);
		TS_ASSERT(cro_data_reader_->map_for_residue_exists(1) );
		TS_ASSERT(!cro_data_reader_->map_for_residue_exists(2) );
		TS_ASSERT_EQUALS(cro_data_reader_->number_of_residues(), 1 );
		std::map <std::string, std::string> expected_map{ {"N1", "N"} };
		TS_ASSERT_EQUALS(cro_data_reader_->get_residue_names_map(1), expected_map);
	}


	void test_access_to_nonexistent_residue() {
		// The map is empty
		try {
			// Mute BACKTRACE message
			set_throw_on_next_assertion_failure();
			// This line should throw the exception if everything works correctly
			cro_data_reader_->get_residue_names_map(1);
			// This line will fail the test if the exception is not thrown
			TS_ASSERT (false);
		}
catch (utility::excn::Exception & e ) {
	// Expected error message
	std::string const error = "The requested residue number is not found.";
	TS_ASSERT_STRING_CONTAINS(e.msg(), error);
}
		// The map exists but there is no information for the requested residue
		std::string meaningful_line("N1   N    1");
		std::istringstream test_stream(meaningful_line);
		cro_data_reader_->create_residue_names_maps(test_stream);
		try {
			// Mute BACKTRACE message
			set_throw_on_next_assertion_failure();
			// This line should throw the exception if everything works correctly
			cro_data_reader_->get_residue_names_map(4);
			// This line will fail the test if the exception is not thrown
			TS_ASSERT (false);
		}
catch (utility::excn::Exception & e ) {
	// Expected error message
	std::string const error = "The requested residue number is not found.";
	TS_ASSERT_STRING_CONTAINS(e.msg(), error);
}

	}

	void test_parse_line_with_wrong_number_of_columns() {
		std::string only_two_columns("N1   N  ");
		std::istringstream test_stream(only_two_columns);
		try {
			// Mute BACKTRACE message
			set_throw_on_next_assertion_failure();
			// This line should throw the exception if everything works correctly
			cro_data_reader_->create_residue_names_maps(test_stream);
			// This line will fail the test if the exception is not thrown
			TS_ASSERT (false);
		}
catch (utility::excn::Exception & e ){
	// Expected error message
	std::string const error = "Wrong number of columns in the parsed line. Please check your input file.";
	TS_ASSERT_STRING_CONTAINS(e.msg(), error);
}
	}

	void test_parse_line_without_number_in_the_third_column() {
		std::string no_number("N1   N  N");
		std::istringstream test_stream(no_number);
		try {
			// Mute BACKTRACE message
			set_throw_on_next_assertion_failure();
			// This line should throw the exception if everything works correctly
			cro_data_reader_->create_residue_names_maps(test_stream);
			// This line will fail the test if the exception is not thrown
			TS_ASSERT (false);
		}
catch (utility::excn::Exception & e ){
	// Expected error message
	std::string const error = "The last column is expected to be a number but it is not. Please check your input file.";
	TS_ASSERT_STRING_CONTAINS(e.msg(), error);
}
	}

	void test_parse_line_with_long_atom_name () {
		std::string long_atom_name("NA12B   N  3");
		std::istringstream test_stream(long_atom_name);
		try {
			// Mute BACKTRACE message
			set_throw_on_next_assertion_failure();
			// This line should throw the exception if everything works correctly
			cro_data_reader_->create_residue_names_maps(test_stream);
			// This line will fail the test if the exception is not thrown
			TS_ASSERT (false);
		}
catch (utility::excn::Exception & e ){
	// Expected error message
	std::string const error = "Atom names are expected to be no longer than 4 characters:";
	TS_ASSERT_STRING_CONTAINS(e.msg(), error);
}
	}

	void test_parse_line_with_four_letter_atom_name () {
		// Four letter atom name should not raise any errors
		std::string four_letter_atom_name("NA12  N  3");
		std::istringstream test_stream(four_letter_atom_name);
		cro_data_reader_->create_residue_names_maps(test_stream);
		TS_ASSERT(cro_data_reader_->map_for_residue_exists(3) );
	}

	void test_nonexistent_file () {
		try {
			// Mute BACKTRACE message
			set_throw_on_next_assertion_failure();
			// This line should throw the exception if everything works correctly
			cro_data_reader_->initialize_from_file("nonexistent_file.txt");
			// This line will fail the test if the exception is not thrown
			TS_ASSERT (false);
		}
catch (utility::excn::Exception & e ){
	// Expected error message
	std::string const error = "Cannot open file ";
	TS_ASSERT_STRING_CONTAINS(e.msg(), error);
}
	}

	void test_initialize_from_file () {
		cro_data_reader_->initialize_from_file("protocols/chromophore/sample_CRO_data_file.txt");
		TS_ASSERT_EQUALS(cro_data_reader_->number_of_residues(), 3 );
		TS_ASSERT_EQUALS(cro_data_reader_->get_residue_names_map(1)["N1"], "N");
		TS_ASSERT_EQUALS(cro_data_reader_->get_residue_names_map(1)["OG1"], "OG1");
		TS_ASSERT_EQUALS(cro_data_reader_->get_residue_names_map(1)["C1"], "C");
		TS_ASSERT_EQUALS(cro_data_reader_->get_residue_names_map(2)["CE1"], "CE1");
		TS_ASSERT_EQUALS(cro_data_reader_->get_residue_names_map(3)["O3"], "O");
	}


	void test_initialize_second_time () {
		// The first initialization
		cro_data_reader_->initialize_from_file("protocols/chromophore/sample_CRO_data_file.txt");
		TS_ASSERT_EQUALS(cro_data_reader_->number_of_residues(), 3 );
		// Attempt to initialize the second time
		std::string meaningful_line("N1   N    1");
		std::istringstream test_stream(meaningful_line);
		try {
			// Mute BACKTRACE message
			set_throw_on_next_assertion_failure();
			// This line should throw the exception if everything works correctly
			cro_data_reader_->create_residue_names_maps(test_stream);
			// This line will fail the test if the exception is not thrown
			TS_ASSERT (false);
		}
catch (utility::excn::Exception & e ){
	// Expected error message
	std::string const error = "One file was already loaded by this object";
	TS_ASSERT_STRING_CONTAINS(e.msg(), error);
}
	}

private:

	protocols::chromophore::ChromophoreDataReaderOP cro_data_reader_;

};
