// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/ysrinivasan/mucintypeglycosylation.cc
/// @brief A protocol to dock glycoprotein glycans in glycosyltransferases at the active site.
/// @author Yashes Srinivasan (yashess@gmail.com)

// Project headers
#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/glycopeptide_docking/GlycopeptideDockingProtocol.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// C++ headers
#include <string>
#include <iostream>

int const SUCCESS( 0 );
int const FAILURE( -1 );

// Main ////////////////////////////////////////////////////////////////////////
int
main( int argc, char *argv[] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		// Initialize Rosetta.
		std::cout << std::endl << "Initializing Rosetta..." << std::endl;
		devel::init( argc, argv );

		// Construct the protocol.
		protocols::glycopeptide_docking::GlycopeptideDockingProtocolOP protocol( new protocols::glycopeptide_docking::GlycopeptideDockingProtocol );

		if ( option[ in::file::native ].active() ) {
			protocol->set_ref_pose_from_filename( option[ in::file::native ] );
		}

		// Distribute the mover.
		protocols::jd2::JobDistributor::get_instance()->go( protocol );

	} catch (utility::excn::Exception const & e ) {
		std::cerr << "Caught exception: " << e.msg() << std::endl;
		return FAILURE;
	}

	return SUCCESS;
}
