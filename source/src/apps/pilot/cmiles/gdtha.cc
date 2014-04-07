// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/cmiles/gdtha.cc
/// @author Christopher Miles (cmiles@uw.edu)

// C/C++ headers
#include <iostream>
#include <map>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>

// Project headers
#include <core/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>

int main(int argc, char* argv[]) {
	try {

  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace core::pose;
  using namespace std;
  devel::init(argc, argv);

  PoseOP ref = core::import_pose::pose_from_pdb(option[OptionKeys::in::file::native]());
  utility::vector1<PoseOP> models = core::import_pose::poseOPs_from_pdbs(option[OptionKeys::in::file::s]());

  map<core::Size, core::Size> all_residues;
  for (core::Size i = 1; i <= ref->total_residue(); ++i) {
    all_residues[i] = i;
  }

  for (utility::vector1<PoseOP>::const_iterator i = models.begin(); i != models.end(); ++i) {
    cout << core::scoring::gdtha(*ref, **i, all_residues) << endl;
  }

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
