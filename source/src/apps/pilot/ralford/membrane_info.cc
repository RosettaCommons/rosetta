// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       apps/pilot/ralford/memrbane_info.cc
///
/// @brief      Unit Test for Membrane Info Object
/// @details
///
/// @author     Rebecca Faye Alford (rfalford12@gmail.com)
/// @note       Last Modified (3/17/14)

// App Headers
#include <devel/init.hh>

// Project Headers
#include <core/membrane/MembraneInfo.hh>
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>

// Package Headers
#include <core/types.hh>

#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using basic::Error;
using basic::Warning;

static basic::Tracer TR( "apps.pilot.ralford.membrane_foldtree" );

/// @brief Load Membrane Pose
core::pose::PoseOP load_pose() {

    using namespace core::import_pose;
    using namespace core::pose;

    TR << "Loading 1afo from PDB" << std::endl;
    PoseOP pose = new Pose();
    pose_from_pdb( *pose, "test/core/membrane/io/1afo_test.pdb" );

    return pose;
}

/// @brief Main Function
int main( int argc, char* argv[] )
{
    try {

		using namespace core::membrane;

        // Initialize Options System, RG, and All Factory_Registrators
        devel::init(argc, argv);

        TR << "Pilot App: Membrane Info Object" << std::endl;
        TR << "Author: Rebecca Alford lm: 3/17/14" << std::endl;
        TR << "Testing membrane info object" << std::endl;

        // Set up a pose from pdb
        core::pose::PoseOP pose = load_pose();

				// Printig pose total residue
				TR << "The number of residues in my pose is " << pose->total_residue() << std::endl;

		// Create embres and membrane data
				utility::vector1< std::pair< int, int > >  embres_map;
        embres_map.resize( 2 );
        embres_map[ 1 ] = std::pair< int, int >( 1, 82 );
        embres_map[ 2 ] = std::pair< int, int >( 41, 83 );

        // Setup the membrane root
        int membrane = 81;

				// Checking that I am passing a valid conformation
				Conformation const & conf = pose->conformation();
				TR << "Printing the size of my conformation " << conf.size() << std::endl;


				TR << "About to dereference" << std::endl;
		// initialize MembraneInfo
		MembraneInfoOP mp = new MembraneInfo( pose->conformation(), embres_map, membrane );

        TR << "Done!" << std::endl;

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cout << "caught exception " << e.msg() << std::endl;
				return -1;
    }
}
