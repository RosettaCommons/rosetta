// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/cmiles/rmsd_partial_thread.cc
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
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>

using core::pose::PDBInfoCOP;
using core::pose::PoseCOP;

int main(int argc, char* argv[]) {
  try {
  using core::Real;
  using core::Size;
  using namespace std;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  devel::init(argc, argv);

  PoseCOP reference = core::import_pose::pose_from_file(option[OptionKeys::in::file::s]()[1], core::import_pose::PDB_file);
  PoseCOP partial_thread = core::import_pose::pose_from_file(option[OptionKeys::in::file::s]()[2], core::import_pose::PDB_file);

  map<Size, Size> residues;

  PDBInfoCOP info = partial_thread->pdb_info();
  for (Size i = 1; i <= partial_thread->size(); ++i) {
    Size reference_pos = info->number(i);
    residues[reference_pos] = i;
  }

  Real rmsd = core::scoring::CA_rmsd(*reference, *partial_thread, residues);
  Real gdtmm = core::scoring::CA_gdtmm(*reference, *partial_thread, residues);
  cout << "rmsd: " << rmsd << ", gdtmm: " << gdtmm << endl;
  } catch ( utility::excn::EXCN_Base const & e ) {
                            std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                }
      return 0;

}
