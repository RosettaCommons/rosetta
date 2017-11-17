// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers
//#include <basic/options/option.hh>

// Utility headers
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>

int main(int argc, char *argv[]) {

	try {

	devel::init(argc,argv);
	core::pose::Pose p,q;
	core::import_pose::pose_from_file(p,basic::options::option[basic::options::OptionKeys::in::file::s]()[1], core::import_pose::PDB_file);
	core::import_pose::pose_from_file(q,basic::options::option[basic::options::OptionKeys::in::file::s]()[2], core::import_pose::PDB_file);
	// core::pose::symmetry::make_symmetric_pose(p);
	// core::pose::symmetry::make_symmetric_pose(q);
	core::Real x = core::scoring::CA_rmsd(p,q);
	std::cout << x << std::endl;
	return 0;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
