// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilat/andre/SymDock.cc
/// @brief  Symmetric Docking protocol

// libRosetta headers
#include <protocols/symmetric_docking/SymDockProtocol.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/moves/MoverContainer.hh>


//#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <utility/excn/Exceptions.hh>

//#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <devel/init.hh>

// C++ headers
//#include <cstdlib>
#include <iostream>
#include <string>

#include <basic/Tracer.hh>

// option key includes

#include <protocols/viewer/viewers.hh>
//#include <basic/options/option_macros.hh>

// add options
//OPT_1GRP_KEY( Boolean, r_SymDock_app, viewer )
// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>


//local options

/*namespace basic{ namespace options{ namespace OptionKeys{
		basic::options::BooleanOptionKey const r_SymDock_app_viewer("enable viewer");
}}}//basic::options::OptionKeys
*/

using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;

void *
SymDock_main_local( void *)
{
	using namespace protocols::symmetric_docking;
	using namespace protocols::simple_moves::symmetry;
	using namespace protocols::jd2;

	SetupForSymmetryMoverOP setup_mover = new SetupForSymmetryMover;
	SymDockProtocolOP dock_mover = new SymDockProtocol;
	protocols::moves::SequenceMoverOP seq_mover = new protocols::moves::SequenceMover;
	seq_mover->add_mover( setup_mover );
	seq_mover->add_mover( dock_mover );
	protocols::jd2::JobDistributor::get_instance()->go( seq_mover );
	return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////


int
main( int argc, char * argv [] )
{
	try{
	using namespace basic::options;
  using namespace basic::options::OptionKeys;
	//	NEW_OPT( r_SymDock_app::viewer, "enable viewer", false );
	//  option.add(r_SymDock_app_viewer, "Enable viewer").def(false);

	devel::init( argc, argv );
	// use viewer if flag given
	if ( option[ OptionKeys::docking::view ]() ){
		protocols::viewer::viewer_main( SymDock_main_local );}
	else{
		SymDock_main_local(NULL);}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
