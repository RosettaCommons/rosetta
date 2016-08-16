// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/polynomial.cxxtest.hh
/// @brief  Test hbond polynomial classes
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/core/init_util.hh>

// Package Headers
#include <core/types.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>

// Project Headers
#include <core/scoring/hbonds/polynomial.hh>
#include <basic/database/open.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

//Auto Headers


class HBPolyTest : public CxxTest::TestSuite {

public:
  void setUp() {
    core_init();
  }

  void tearDown(){}


  void test_Polynomial_1d_model(){

    test::UTracer UT("core/scoring/hbonds/Polynomial_1d_model.u");

    UT << "let p(x) = 5x^2 - 1" << std::endl;
    utility::vector1< core::Real > coefs;
    coefs.push_back(5);
		coefs.push_back(0);
    coefs.push_back(-1);

    core::scoring::hbonds::Polynomial_1d p(
      "test_polynomial",
      core::scoring::hbonds::hbgd_AHdist,
      -3 /*xmin*/,
      3 /*xmax*/,
      10.0 /*min_val*/,
      10.0 /*max_val*/,
      -.2 /*root1*/,
      .2 /*root2*/,
      3 /*degree*/,
      coefs);

		core::Real value, deriv;
    p(0, value, deriv);
		UT << "p(0) = " << value << ", p'(0) = " << deriv << std::endl;
    TS_ASSERT(value == -1);
    TS_ASSERT(deriv == 0);

    p(2, value, deriv);
		UT << "p(2) = " << value << ", p'(2) = " << deriv << std::endl;
    TS_ASSERT(value == 19);
    TS_ASSERT(deriv == 20);

    UT << p;
  }

	void test_polynomial_1d_check_invariants(){

		utility::vector1< core::Real > coefs;
		coefs.push_back(5);
		coefs.push_back(0);
		coefs.push_back(-1);

		try{
			// xmax < xmin
			core::scoring::hbonds::Polynomial_1d p(
				"test_polynomial",
				core::scoring::hbonds::hbgd_AHdist,
				-3 /*xmin*/,
				-4 /*xmax*/,
				10.0 /*min_val*/,
				10.0 /*max_val*/,
				-.2 /*root1*/,
				.2 /*root2*/,
				3 /*degree*/,
				coefs);
			TS_ASSERT(false);
		} catch {}

		try{
			// max_val < min_val
			core::scoring::hbonds::Polynomial_1d p(
				"test_polynomial",
				core::scoring::hbonds::hbgd_AHdist,
				-3 /*xmin*/,
				3 /*xmax*/,
				10.0 /*min_val*/,
				-10.0 /*max_val*/,
				-.2 /*root1*/,
				.2 /*root2*/,
				3 /*degree*/,
				coefs);
			TS_ASSERT(false);
		} catch {}

		try{
			// coefficients are empty
			utility::vector1< core::Real > empty_coefs;

			core::scoring::hbonds::Polynomial_1d p(
				"test_polynomial",
				core::scoring::hbonds::hbgd_AHdist,
				-3 /*xmin*/,
				3 /*xmax*/,
				10.0 /*min_val*/,
				-10.0 /*max_val*/,
				-.2 /*root1*/,
				.2 /*root2*/,
				3 /*degree*/,
				empty_coefs);
			TS_ASSERT(false);
		} catch {}

		try{
			// degree != degree of passed in coefficients
			core::scoring::hbonds::Polynomial_1d p(
				"test_polynomial",
				core::scoring::hbonds::hbgd_AHdist,
				-3 /*xmin*/,
				3 /*xmax*/,
				10.0 /*min_val*/,
				-10.0 /*max_val*/,
				-.2 /*root1*/,
				.2 /*root2*/,
				2 /*degree*/,
				empty_coefs);
			TS_ASSERT(false);
		} catch {}

	}

	void test_read_polynomials_from_database(){

		test::UTracer UT("core/scoring/hbonds/read_polynomials_from_database.u");

		UT << "Read all the polynomials from the hbond database" << std::endl;
		core::scoring::hbonds::HBondOptionsCOP hb_options( new core::scoring::hbonds::HBondOptions() );
		core::scoring::hbonds::HBondDatabaseCOP hb_database( core::scoring::hbonds::HBondDatabase::get_database(hb_options->params_database_tag()));

		std::stringstream HBPoly1D_fname;
		HBPoly1D_fname  << "scoring/score_functions/hbonds/" << hb_options->params_database_tag() << "/HBPoly1D.csv";
		utility::io::izstream s;
		basic::database::open(s, HBPoly1D_fname.str());
		std::string line;
		core::scoring::hbonds::Polynomial_1dCOP p;
		while (getline( s, line)){
			std::istringstream l(line);
			utility::vector1<std::string> tokens;
			tokens = utility::string_split( line, ',');
			std::string polynomial_name(tokens[2]);
			UT << *(hb_database->HBPoly1D_from_name(polynomial_name)) << std::endl;
			UT << hb_database->HBPoly1D_from_name(polynomial_name)->show_values() << std::endl;
		}
	}

};
