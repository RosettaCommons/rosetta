// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/cmiles/dssp.cc
/// @author Christopher Miles (cmiles@uw.edu)

// C/C++ headers
#include <iostream>

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
#include <core/scoring/dssp/Dssp.hh>

int main(int argc, char* argv[]) {
	try {

  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  devel::init(argc, argv);
  core::pose::PoseCOP pose = core::import_pose::pose_from_file(option[OptionKeys::in::file::s]()[1], core::import_pose::PDB_file);
  core::pose::PDBInfoCOP info = pose->pdb_info();

  core::scoring::dssp::Dssp dssp(*pose);
  dssp.dssp_reduced();

  std::cout << "Rosetta_Residue" << "\t" << "PDB_Residue" << "\t" << "Secondary_Structure" << std::endl;
  for (core::Size i = 1; i <= pose->total_residue(); ++i) {
    std::cout << i << "\t" << info->number(i) << "\t" << dssp.get_dssp_secstruct(i) << std::endl;
  }

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
