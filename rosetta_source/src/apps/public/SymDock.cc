// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file /src/apps/pilat/andre/SymDock.cc
/// @brief  Symmetric Docking protocol

// libRosetta headers
#include <protocols/symmetric_docking/SymDockProtocol.hh>
#include <protocols/moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/moves/MoverContainer.hh>


//#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>

#include <protocols/jd2/JobDistributor.hh>

//#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <devel/init.hh>

// C++ headers
//#include <cstdlib>
#include <iostream>
#include <string>

#include <basic/Tracer.hh>

// option key includes

using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;

int
main( int argc, char * argv [] )
{
	using namespace protocols::symmetric_docking;
	using namespace protocols::moves::symmetry;
	using namespace protocols::jd2;

	// initialize core
	devel::init(argc, argv);

	protocols::symmetric_docking::SymDock_main();
}
