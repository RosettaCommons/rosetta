// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/pilot/jianqing/ncaa_crosslink.cc
/// @brief
/// @NCAA utility
/// @author Jianqing Xu (xubest@gmail.com)
/// 09/09/2011


#include <protocols/ncaa_crosslink/NcaaPreCrosslink.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <devel/init.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>
//#include <utility/tools/make_vector1.hh>

// option key includes
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/option.hh>


int
main( int argc, char * argv [] )
{

	try {

	using namespace basic::options;
	using namespace protocols::ncaa_crosslink;
	using namespace protocols::jd2;

//	NcaaPreCrosslink::register_options();
	// initialize core
//	devel::init(argc, argv);

//	NcaaPreCrosslinkOP npc;
//	npc = new NcaaPreCrosslink();

 	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}


}


/*

int main( int argc, char * argv [] ) {
	devel::init(argc, argv);

	core::pose::PoseOP pose = new core::pose::Pose();
	core::import_pose::pose_from_pdb( *pose, "./5A_low_1");


	core::scoring::ScoreFunctionCOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "mm_std") ;


	core::pack::task::TaskFactoryOP tf = new core::pack::task::TaskFactory;

	tf->push_back( new core::pack::task::operation::RestrictToRepacking );
	tf->push_back( new core::pack::task::operation::InitializeFromCommandline );
	tf->push_back( new core::pack::task::operation::IncludeCurrent );
	tf->push_back( new core::pack::task::operation::NoRepackDisulfides );


	utility::vector1< core::Size > const movable_jumps(1);

//	protocols::docking::DockJumps const movable_jumps  =  utility::vector1<core::Size> ;
//	tf->push_back( new RestrictToInterface( movable_jumps() ) );


//	protocols::docking::setup_foldtree( pose, "LH_J", movable_jumps );
}


*/

