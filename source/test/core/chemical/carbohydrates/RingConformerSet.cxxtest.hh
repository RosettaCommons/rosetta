// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	 RingConformerSet.cxxtest.hh
/// @brief   Test suite for ring conformer set building and associated methods
/// @author  Labonte

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/chemical/carbohydrates/RingConformerSet.hh>

// Basic headers
//#include <basic/options/option.hh>
//#include <basic/options/keys/out.OptionKeys.gen.hh>

// Utility header
#include <utility/vector1.hh>


class RingConformerSetTests : public CxxTest::TestSuite {
public:
	// Standard methods ///////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		//using namespace basic::options;
		using namespace core::chemical::carbohydrates;

		core_init_with_additional_options("-out:levels core.chemical.carbohydrates.RingConformerSet:400");

		//option[OptionKeys::out::levels](StringVectorOptionKey("core.chemical.carbohydrates.RingConformerSet:400"));

		set5_ = new RingConformerSet(5);
		set6_ = new RingConformerSet(6);
	}

	// Destruction
	void tearDown()
	{}

	// Tests //////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Confirm that the proper number of conformers are loaded into 5- and 6-membered ring sets.
	void test_get_all_nondegenerate_conformers()
	{
		TS_TRACE("Testing get_all_nondegenerate_conformers() method of RingConformerSet for 5- and 6-membered rings.");
		TS_ASSERT_EQUALS(set5_->get_all_nondegenerate_conformers().size(), 20);
		TS_ASSERT_EQUALS(set6_->size(), 38);
	}

	// Confirm that the proper conformers are loaded by get_ideal_conformer_by_name() and
	// get_ideal_conformer_by_CP_parameters.
	void test_get_ideal_conformer_by_descriptor_methods()
	{
		using namespace core;
		using namespace core::chemical::carbohydrates;
		using namespace utility;

		TS_TRACE("Testing get_ideal_conformer_by_name() and get_ideal_conformer_by_CP_parameters() methods of"
				" RingConformerSet for 5- and 6-membered rings.");

		vector1<Real> params5, params6, params_bad;
		params5.resize(2);
		params6.resize(3);

		params5[q] = 0.4;
		params5[PHI] = 180.0;

		params6[q] = 0.55;
		params6[PHI] = 180.0;
		params6[THETA] = 90.0;

		TS_ASSERT_EQUALS(set5_->get_ideal_conformer_by_name("EO"),
				set5_->get_ideal_conformer_by_CP_parameters(params5));
		TS_ASSERT_EQUALS(set6_->get_ideal_conformer_by_name("BO,3"),
				set6_->get_ideal_conformer_by_CP_parameters(params6));

		// Test that rounding is handled properly.
		params5[PHI] = 43.2;  // should round to 36.0

		TS_ASSERT_EQUALS(set5_->get_ideal_conformer_by_name("E1"),
				set5_->get_ideal_conformer_by_CP_parameters(params5));

		params6[PHI] = 123.4;  // should round to 120.0
		params6[THETA] = 123.4;  // should round to 135.0

		TS_ASSERT_EQUALS(set6_->get_ideal_conformer_by_name("5E"),
				set6_->get_ideal_conformer_by_CP_parameters(params6));

		// Test that chairs are handled correctly, as they reside at the poles where phi is meaningless.
		params6[THETA] = 0.0;

		TS_ASSERT_EQUALS(set6_->get_ideal_conformer_by_name("4C1"),
				set6_->get_ideal_conformer_by_CP_parameters(params6));

		// Test for bad input.
		TS_ASSERT_EQUALS(set5_->get_ideal_conformer_by_name("FOO"), NULL);
		TS_ASSERT_EQUALS(set6_->get_ideal_conformer_by_CP_parameters(params_bad), NULL);
		params_bad.push_back(0.0);  // still bad because not enough params for a 6-membered ring
		TS_ASSERT_EQUALS(set6_->get_ideal_conformer_by_CP_parameters(params_bad), NULL);
		params_bad.push_back(180.0);
		params_bad.push_back(90.0);  // still bad because q = 0.0 (planar)
		TS_ASSERT_EQUALS(set6_->get_ideal_conformer_by_CP_parameters(params_bad), NULL);
		params_bad.push_back(90.0);  // still bad because too many params
		TS_ASSERT_EQUALS(set6_->get_ideal_conformer_by_CP_parameters(params_bad), NULL);
	}

private:
	// Private data ////////////////////////////////////////////////////////////
	core::chemical::carbohydrates::RingConformerSetOP set5_;
	core::chemical::carbohydrates::RingConformerSetOP set6_;

};  // class RingConformerSetTests
