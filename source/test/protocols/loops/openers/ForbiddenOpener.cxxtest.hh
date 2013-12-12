// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Unit test for the ForbiddenOpener class.
/// @author Kale Kundert (kalekundert@ucsf.edu)

// Test headers
#include <cxxtest/TestSuite.h>

// Project headers
#include <protocols/loops/openers/Opener.hh>
#include <protocols/loops/openers/ForbiddenOpener.hh>

using namespace protocols::loops;

class ForbiddenOpenerTests : public CxxTest::TestSuite {

 public:
  void test_constructor() {
		openers::OpenerOP opener = new openers::ForbiddenOpener();
  }

};

