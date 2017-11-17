// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/pilat/andre/SymDock.cc
/// @brief  Symmetric Docking protocol

// libRosetta headers
#include <protocols/symmetric_docking/SymDockProtocol.hh>

#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>


// option key includes

using basic::Error;
using basic::Warning;

using namespace core;

int
main( int argc, char * argv [] )
{
	try {

		using namespace protocols::symmetric_docking;

		// initialize core
		devel::init(argc, argv);

		protocols::symmetric_docking::SymDock_main();

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
