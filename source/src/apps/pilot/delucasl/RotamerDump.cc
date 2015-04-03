// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/public/design/fixbb.cc
/// @brief  Fixed backbone design.  Can do side chain minimization after PackRotamers by using the flag -minimize_sidechains.  This is SLOW.

//core library
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

#include <basic/options/option.hh>

#include <protocols/RotamerDump/RotamerDumpMover.hh>

//utilities
#include <protocols/jd2/JobDistributor.hh>
#include <devel/init.hh>

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <utility/vector1.hh>

//Auto Headers


//local options
namespace basic{ namespace options{ namespace OptionKeys{
basic::options::BooleanOptionKey const minimize_sidechains("minimize_sidechains");
}}}//basic::options::OptionKeys

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
    try {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	option.add( minimize_sidechains, "Do minimization of side chains after rotamer packing").def(false);

	devel::init(argc, argv);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// end of setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//create a task factory: this will create a new PackerTask for each input pose
	core::pack::task::TaskFactoryOP main_task_factory = new core::pack::task::TaskFactory;
	main_task_factory->push_back( new core::pack::task::operation::InitializeFromCommandline );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory->push_back( new core::pack::task::operation::ReadResfile );
	}

	//create a ScoreFunction from commandline options
	core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();

	protocols::moves::MoverOP rotamer_mover(new protocols::RotamerDump::RotamerDumpMover(main_task_factory,score_fxn));


	protocols::jd2::JobDistributor::get_instance()->go(rotamer_mover);

    } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1
	}
	return 0;
}
