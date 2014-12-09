// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/OptionCollection.cxxtest.hh
/// @brief  test suite for utility::OptionCollection
/// @author Matt O'Meara (mattjomeara@gmail.com)

// Package headers
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/IntegerOptionKey.hh>

class OptionCollectionTests : public CxxTest::TestSuite {
public:
	void
	setup() {}


	void test_check_specs() {
		utility::options::OptionCollection option;
		utility::options::IntegerOptionKey test_integer_key("id_a", "identifier_a", "code_a");
		option.add( test_integer_key, "test option key" ).def(4);
		option.check_specs();
	}
};


