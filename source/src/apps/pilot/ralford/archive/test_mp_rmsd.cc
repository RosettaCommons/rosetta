// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/pilot/ralford/test_mp_rmsd.cc
///
/// @brief  Trial modified centroid rigid body cycles
///    Last Modified: 10/22/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/membrane/metrics.hh>

#include <core/import_pose/import_pose.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>

using namespace protocols::membrane;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.ralford.membrane_rmsd" );

/// @brief Main method
int
main( int argc, char * argv [] )
{
	using namespace protocols::membrane;
	using namespace core::import_pose;

	try {

		// Devel init factories
		devel::init(argc, argv);

		TR << "Hello to loading in poses" << std::endl;

		core::pose::PoseOP native_pose = new Pose();
		pose_from_file( *native_pose, "1afo_in.pdb" , core::import_pose::PDB_file);

		core::pose::PoseOP test_pose = new Pose();
		pose_from_file( *test_pose, "1afo_decoy.pdb" , core::import_pose::PDB_file);

		// Calculate rmsd
		//TR << "RMSD is: " << compute_mp_rmsd_with_super( *native_pose, *test_pose ) << std::endl;


		return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
