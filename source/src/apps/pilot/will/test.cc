// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
//#include <basic/options/option.hh>

// Utility headers
#include <core/init/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>

int main(int argc, char *argv[]) {	

	try {

	core::init::init(argc,argv);
	core::pose::Pose p,q;
	core::import_pose::pose_from_pdb(p,basic::options::option[basic::options::OptionKeys::in::file::s]()[1]);
	core::import_pose::pose_from_pdb(q,basic::options::option[basic::options::OptionKeys::in::file::s]()[2]);
	// core::pose::symmetry::make_symmetric_pose(p);
	// core::pose::symmetry::make_symmetric_pose(q);	
	core::Real x = core::scoring::CA_rmsd(p,q);
	std::cout << x << std::endl;
	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}
