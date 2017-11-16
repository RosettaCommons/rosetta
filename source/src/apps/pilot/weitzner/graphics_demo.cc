// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/weitzner/long_loop.cc
/// @brief
/// @author Brian Weitzner (brian.weitzner@jhu.edu)

// Unit Headers

// Package Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/CCDLoopClosureMover.hh>

#include <protocols/moves/PyMOLMover.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <protocols/forge/methods/util.hh>

// Numeric Headers and ObjexxFCL Headers

// C++ headers
#include <string>

#include <protocols/viewer/viewers.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

//#include <protocols/viewer/ConformationViewer.hh>

using namespace core;
using namespace protocols::loops;
using namespace protocols::forge::methods;

static basic::Tracer TR( "pilot_apps.weitzner.long_loop" );

LoopsOP
read_loops_from_text_file() {
	LoopsOP loops = new Loops();

	//The following will be replaced by code to read in loop definitions from a text file
	Size loop1_start = 98;
	Size loop1_end = 112;
	Size loop1_cutpoint = 105;

	// There will be a for loop to continually add the loops to the loop set
	loops->add_loop( Loop(loop1_start, loop1_end, loop1_cutpoint) );
	return loops;
}

void
set_all_loop_dihedrals_to_180( pose::PoseOP pose, const Loops::LoopList & loop_list) {
	TR << "Setting all loop dihedrals to 180 degrees (opening the loop)." << std::endl;

	Real my_dihedral = 180.;
	for ( Size i=loop_list[1].start(); i <= loop_list[1].stop(); i++ ) {
		pose->set_phi(i, my_dihedral);
		pose->set_psi(i, my_dihedral);
	}
}

void*
my_main( void* ) {
	pose::PoseOP pose = new pose::Pose();

	core::import_pose::pose_from_file(*pose, "src/apps/pilot/weitzner/1bzq.pdb");// changing to relative path, so it could work from diff places... "/work/rosetta/rosetta_source/src/apps/pilot/weitzner/1bzq.pdb", core::import_pose::PDB_file);

	TR << "The pose's sequence is: " << pose->sequence() << std::endl;

	/*
	protocols::viewer::ConformationViewerOP cv = new protocols::viewer::ConformationViewer();

	cv->attach_to(pose->conformation());
	cv->display_func();
	*/
	protocols::viewer::add_conformation_viewer( pose->conformation(), "1bzq.pdb", 450, 450, true );

	// Hey... let's check how its actually looks like... →→→ ☆★PyMOL★☆
	protocols::moves::AddPyMOLLink(*pose, true); // Lets ask PyMOL to store history...

	// Eventually this will take a text file as input
	LoopsOP loops = read_loops_from_text_file();

	pose->fold_tree( fold_tree_from_loops( *pose, *loops ) );
	TR << pose->fold_tree();
	set_all_loop_dihedrals_to_180( pose, loops->loops() );

	kinematics::MoveMapOP mm = new kinematics::MoveMap();
	mm->set_bb(false);
	loops->loops()[1].switch_movemap( *mm, id::BB, true);

	CCDLoopClosureMover ccd  = CCDLoopClosureMover( loops->loops()[1], mm);
	ccd.apply(*pose);
	pose->dump_pdb("src/apps/pilot/weitzner/1bzq_closed.pdb");
	return 0;

}

int
main( int argc, char * argv [] ) {

	try {

		devel::init(argc, argv);
		protocols::viewer::viewer_main( my_main );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


