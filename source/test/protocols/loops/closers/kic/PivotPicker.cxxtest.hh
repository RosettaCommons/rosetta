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
#include <protocols/loops/Loop.hh>
#include <protocols/loops/closers/kic/PivotPicker.hh>

using namespace protocols::loops::Loop;
using namespace protocols::loops::closers::kic::PivotPicker;

class PivotPickerTests : public CxxTest::TestSuite {

 public:

  void test_check_pivot_residues() {
		Loop loop = Loop(2, 5);
		PivotPicker picker;

		picker.check_pivot_residues(2, 3, 5, loop);  // Ok.
		picker.check_pivot_residues(1, 3, 5, loop);  // Not ok.
		picker.check_pivot_residues(2, 3, 6, loop);  // Not ok.
		picker.check_pivot_residues(2, 5, 3, loop);  // Not ok.
  }

  void test_check_pivot_atoms() {
		Loop loop = Loop(2, 5);
		PivotPicker picker;

		picker.check_pivot_atoms(4, 9, 15);  // Ok.
		picker.check_pivot_atoms(3, 9, 15);  // Not ok.
		picker.check_pivot_atoms(4, 9, 16);  // Not ok.
		picker.check_pivot_atoms(4, 15, 9);  // Not ok.
  }

};

