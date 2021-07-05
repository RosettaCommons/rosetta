// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/chromophore/ReadResidueCoordinatesFromPDBTests.cxxtest.hh
/// @brief  unit tests for ReadResidueCoordinatesFromPDB class
/// @author Nina Bozhanova (nbozhanova@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/chromophore/ReadResidueCoordinatesFromPDB.hh>

// Core Headers

// Utility, etc Headers
#include <basic/Tracer.hh>

#include <numeric/xyzVector.hh> // AUTO IWYU For xyzVector

static basic::Tracer TR("ReadResidueCoordinatesFromPDBTests");


class ReadResidueCoordinatesFromPDBTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		res_coordinates_reader_ = protocols::chromophore::ReadResidueCoordinatesFromPDBOP (new protocols::chromophore::ReadResidueCoordinatesFromPDB());
	}

	void tearDown(){
		res_coordinates_reader_ = 0;
	}

	void test_no_coordinates() {
		TS_ASSERT_EQUALS(res_coordinates_reader_->number_of_residues (), 0);
		// Try to get coordinates when they do not exist
		try {
			// Mute BACKTRACE message
			set_throw_on_next_assertion_failure();
			// This line should throw the exception if everything works correctly
			res_coordinates_reader_->get_residue_coordinates(6, 'A');
			// This line will fail the test if the exception is not thrown
			TS_ASSERT (false);
		}
catch (utility::excn::Exception & e ){
	// Expected error message
	std::string const error = "No coordinates are available for the residue.";
	TS_ASSERT_STRING_CONTAINS(e.msg(), error);
}

	}

	void test_read_from_nonexistent_file() {
		try {
			// Mute BACKTRACE message
			set_throw_on_next_assertion_failure();
			// This line should throw the exception if everything works correctly
			utility::vector1 <std::tuple <int, char> > residues;
			residues.push_back(std::tuple <int, char>(64, 'A'));
			res_coordinates_reader_->read_coordinates_from_file("nonexistent_file.txt", residues);
			// This line will fail the test if the exception is not thrown
			TS_ASSERT (false);
		}
catch (utility::excn::Exception & e ){
	// Expected error message
	std::string const error = "Cannot open file ";
	TS_ASSERT_STRING_CONTAINS(e.msg(), error);
}
	}

	void test_read_from_file() {
		// Initializing options system
		core_init();

		utility::vector1 <std::tuple <int, char> > residues;
		residues.push_back(std::tuple <int, char>(34, 'A'));
		residues.push_back(std::tuple <int, char>(35, 'A'));
		res_coordinates_reader_->read_coordinates_from_file("protocols/chromophore/test.pdb", residues);
		// Was reading of the coordinates successful?
		TS_ASSERT(res_coordinates_reader_->coordinates_exist(34, 'A'));
		TS_ASSERT(res_coordinates_reader_->coordinates_exist(35, 'A'));
		TS_ASSERT_EQUALS(res_coordinates_reader_->number_of_residues (), 2);
		// Expected coordinates for residue #1
		std::string res1_expected_atom_name_first = "N";
		double res1_x_coordinate_first = 7.132;
		double res1_y_coordinate_first = -8.342;
		double res1_z_coordinate_first = 23.730;
		std::string res1_expected_atom_name_last = "2HG";
		double res1_x_coordinate_last = 8.520;
		double res1_y_coordinate_last = -9.980;
		double res1_z_coordinate_last = 25.336;
		// Expected coordinates for residue #2
		std::string res2_expected_atom_name = "C";
		double res2_x_coordinate = 9.278;
		double res2_y_coordinate = -3.012;
		double res2_z_coordinate = 25.323;
		// Returned coordinates
		utility::vector1 <std::tuple <std::string, core::Vector> > res1_data = res_coordinates_reader_->get_residue_coordinates(34, 'A');
		utility::vector1 <std::tuple <std::string, core::Vector> > res2_data = res_coordinates_reader_->get_residue_coordinates(35, 'A');
		// Checking residue #1
		TS_ASSERT_EQUALS(res1_data.size(), 15);
		TS_ASSERT_EQUALS(std::get<0>(res1_data[1]), res1_expected_atom_name_first);
		TS_ASSERT_DELTA(std::get<1>(res1_data[1]).x(), res1_x_coordinate_first, 0.0001);
		TS_ASSERT_DELTA(std::get<1>(res1_data[1]).y(), res1_y_coordinate_first, 0.0001);
		TS_ASSERT_DELTA(std::get<1>(res1_data[1]).z(), res1_z_coordinate_first, 0.0001);
		TS_ASSERT_EQUALS(std::get<0>(res1_data.back()), res1_expected_atom_name_last);
		TS_ASSERT_DELTA(std::get<1>(res1_data.back()).x(), res1_x_coordinate_last, 0.0001);
		TS_ASSERT_DELTA(std::get<1>(res1_data.back()).y(), res1_y_coordinate_last, 0.0001);
		TS_ASSERT_DELTA(std::get<1>(res1_data.back()).z(), res1_z_coordinate_last, 0.0001);
		// Checking residue #2
		TS_ASSERT_EQUALS(res2_data.size(), 7);
		TS_ASSERT_EQUALS(std::get<0>(res2_data[3]), res2_expected_atom_name);
		TS_ASSERT_DELTA(std::get<1>(res2_data[3]).x(), res2_x_coordinate, 0.0001);
		TS_ASSERT_DELTA(std::get<1>(res2_data[3]).y(), res2_y_coordinate, 0.0001);
		TS_ASSERT_DELTA(std::get<1>(res2_data[3]).z(), res2_z_coordinate, 0.0001);

	}

	void test_read_from_file_nonexistent_residue() {
		// Initializing options system
		core_init();
		utility::vector1 <std::tuple <int, char> > residues;
		residues.push_back(std::tuple <int, char>(34, 'B'));
		try {
			// Mute BACKTRACE message
			set_throw_on_next_assertion_failure();
			// This line should throw the exception if everything works correctly
			res_coordinates_reader_->read_coordinates_from_file("protocols/chromophore/test.pdb", residues);
			// This line will fail the test if the exception is not thrown
			TS_ASSERT (false);
		}
catch (utility::excn::Exception & e ){
	// Expected error message
	std::string const error = "No residues fulfill the specified criteria.";
	TS_ASSERT_STRING_CONTAINS(e.msg(), error);
}
	}

	void test_read_from_file_one_residue_do_not_exist() {
		// Initializing options system
		core_init();
		utility::vector1 <std::tuple <int, char> > residues;
		residues.push_back(std::tuple <int, char>(34, 'B'));
		residues.push_back(std::tuple <int, char>(34, 'A'));
		try {
			// Mute BACKTRACE message
			set_throw_on_next_assertion_failure();
			// This line should throw the exception if everything works correctly
			res_coordinates_reader_->read_coordinates_from_file("protocols/chromophore/test.pdb", residues);
			// This line will fail the test if the exception is not thrown
			TS_ASSERT (false);
		}
catch (utility::excn::Exception & e ){
	// Expected error message
	std::string const error = "The number of found residues does not match the number of requested residues.";
	TS_ASSERT_STRING_CONTAINS(e.msg(), error);
}
	}

	void test_read_from_stream() {
		// Initializing options system
		core_init();
		std::string pdb_fragment("ATOM     18 2HD2 LEU A  62       6.321   3.158  -1.896  1.00 14.96           H\n"
			"ATOM     19 3HD2 LEU A  62       6.854   1.706  -1.646  1.00 14.96           H\n"
			"ATOM     20  N   THR A  63       1.843   4.013   1.462  1.00  0.00           N\n"
			"ATOM     21  CA  THR A  63       0.836   4.998   1.855  1.00  0.00           C\n");
		std::istringstream test_stream(pdb_fragment);
		utility::vector1 <std::tuple <int, char> > residues;
		residues.push_back(std::tuple <int, char>(63, 'A'));
		res_coordinates_reader_->read_coordinates (test_stream, residues);
		// Expected coordinates
		std::string expected_atom_name_1 = "N";
		double x_coordinate_1 = 1.843;
		double y_coordinate_1 = 4.013;
		double z_coordinate_1 = 1.462;
		std::string expected_atom_name_2 = "CA";
		double x_coordinate_2 = 0.836;
		double y_coordinate_2 = 4.998;
		double z_coordinate_2 = 1.855;
		// Returned coordinates
		utility::vector1 <std::tuple <std::string, core::Vector> > res_data = res_coordinates_reader_->get_residue_coordinates(63, 'A');
		TS_ASSERT_EQUALS(res_data.size(), 2);
		TS_ASSERT_EQUALS(std::get<0>(res_data[1]), expected_atom_name_1);
		TS_ASSERT_DELTA(std::get<1>(res_data[1]).x(), x_coordinate_1, 0.0001);
		TS_ASSERT_DELTA(std::get<1>(res_data[1]).y(), y_coordinate_1, 0.0001);
		TS_ASSERT_DELTA(std::get<1>(res_data[1]).z(), z_coordinate_1, 0.0001);
		TS_ASSERT_EQUALS(std::get<0>(res_data[2]), expected_atom_name_2);
		TS_ASSERT_DELTA(std::get<1>(res_data[2]).x(), x_coordinate_2, 0.0001);
		TS_ASSERT_DELTA(std::get<1>(res_data[2]).y(), y_coordinate_2, 0.0001);
		TS_ASSERT_DELTA(std::get<1>(res_data[2]).z(), z_coordinate_2, 0.0001);
	}

	void test_try_save_residue_information_twice() {
		// Initializing options system
		core_init();
		utility::vector1 <std::tuple <int, char> > residues;
		residues.push_back(std::tuple <int, char>(33, 'A'));
		// Read coordinates
		res_coordinates_reader_->read_coordinates_from_file("protocols/chromophore/test.pdb", residues);
		// Was reading of the coordinates successful?
		TS_ASSERT(res_coordinates_reader_->coordinates_exist(33, 'A'));
		// Try to read coordinates of the residue the second time
		try {
			// Mute BACKTRACE message
			set_throw_on_next_assertion_failure();
			// This line should throw the exception if everything works correctly
			res_coordinates_reader_->read_coordinates_from_file("protocols/chromophore/test.pdb", residues);
			// This line will fail the test if the exception is not thrown
			TS_ASSERT (false);
		}
catch (utility::excn::Exception & e ){
	// Expected error message
	std::string const error = "The coordinates for this residue seems to be already saved";
	TS_ASSERT_STRING_CONTAINS(e.msg(), error);
}
	}


private:

	protocols::chromophore::ReadResidueCoordinatesFromPDBOP res_coordinates_reader_;


};
