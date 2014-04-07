// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/dan/LoopExtend.cc
/// @brief Loop Extend protocol
/// @author Daniel J. Mandell

// Unit Headers
#include <devel/loop_extend/LoopExtendMover.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>

#include <protocols/loops/Loops.hh>

// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh> //so that the job distributor call can increase the OP count

//for debugging
// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>
// AUTO-REMOVED #include <protocols/loops/loops_main.hh>

// Utility Headers
#include <devel/init.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <protocols/jd2/JobDistributor.hh>
#include <utility/exit.hh>

// Project headers
// Auto-header: duplicate removed #include <protocols/loops/Loops.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>

using basic::T;
using basic::Error;
using basic::Warning;

//replaces cout
static basic::Tracer tr("apps.pilot.dan.LoopExtend");

///@brief main method for loop extension.
int
main( int argc, char* argv[] )
{
    try {
	using protocols::loops::Loops;
	using protocols::loops::Loop;
	devel::init(argc, argv);

	// read pdb file
	core::pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, core::options::option[ basic::options::OptionKeys::loops::input_pdb ]().name() );

	// read loops file
	protocols::loops::Loops loops;
	loops.read_loop_file( basic::options::option[ basic::options::OptionKeys::loops::loop_file ].value()[1] );
	tr << "initial loops " << loops << std::flush;

	//check that cutpoint is not outside the loop and if necessary notify user that only first loop is used
	Loop extend_loop( *(loops.v_begin()) ); // make a local copy of the loop

	if ( (extend_loop.cut() < extend_loop.start()) || (extend_loop.cut() >= extend_loop.stop() ) ) {
		tr << "Problem with loop: " << extend_loop << std::endl;
		utility_exit_with_message("Cutpoint must be between loop start and end-1 for extension");
	}
	/*
	for( Loops::iterator it=loops.v_begin(), it_end=loops.v_end(); it != it_end; ++it ){
		if( (it->cut() < it->start()) || (it->cut() >= it->stop()) ) { // the cutpoint is not within the loop
			Loop loop( *it );
			tr << "Problem with loop: " << loop << std::endl;
			utility_exit_with_message("Cutpoint must be between loop start and end-1 for extension");
		}

		if (it != loops.v_begin()) {
			tr << "Ignoring loop: " << it.begin() << std::endl;
			tr << "Note: LoopExtend protocol only uses first loop definition in loop file" << std::endl;
		}

	}
	*/
	//get the loop to extend
	//protocols::loops::Loop loop=loops.begin();

	// get the length of the insertion
	core::Size extend_len = basic::options::option[ basic::options::OptionKeys::loops::extend_length ]();
	// do the insertion
	devel::loop_extend::LoopExtendMover extend_mover( extend_loop, extend_len );
	extend_mover.apply(pose);
	// get the output filename
	std::string outfile(
						basic::options::option[ basic::options::OptionKeys::out::path::path ]().name() + "/" +
						basic::options::option[ basic::options::OptionKeys::out::prefix ] + ".pdb"
						);
	pose.dump_pdb(outfile);
	tr << "success" << std::endl;

    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;
}
