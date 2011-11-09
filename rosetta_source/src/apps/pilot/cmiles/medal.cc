// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/cmiles/medal.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Utility headers
#include <devel/init.hh>

// Project headers
#include <protocols/medal/MedalMain.hh>
#include <protocols/viewer/viewers.hh>

int main(int argc, char* argv[]) {
  devel::init(argc, argv);
  protocols::viewer::viewer_main(
      protocols::medal::Medal_main);
}
