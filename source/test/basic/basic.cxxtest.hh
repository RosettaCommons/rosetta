// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/basic.cxxtest.hh
/// @brief  basic.cxxtest: placeholder for main() function when creating unit test executable
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifndef INCLUDED_basic_basic_cxxtest_HH
#define INCLUDED_basic_basic_cxxtest_HH

// When compiling tests for CxxTest using the --part argument to the Python script, no main()
// method is added to the generated .cpp files. But there needs to be one file that has a main.
// Let's assume inside each "project" there will be a project.cxxtest.hh which has no tests or anything
// in it. This file will get processed with the --root switch which adds a main() method to the
// generated .cpp file.  This is that file for the basic project.


#endif INCLUDED_basic_cxxtest_HH


