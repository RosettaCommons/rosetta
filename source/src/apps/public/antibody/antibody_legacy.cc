// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file Aroop Sircar ( aroopsircar@yahoo.com )
/// @brief Legacy RosettaAntibody application as written by Arvind Sivasubramanian and Aroop Sircar, and ported to Rosetta3


#include <protocols/jobdist/standard_mains.hh>

// Rosetta Headers
#include <basic/Tracer.hh>
#include <devel/init.hh>

#include <protocols/antibody_legacy/AntibodyModeler.hh>
#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


////////////////////////////////////////////////////////
using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;

int
main( int argc, char * argv [] )
{
	try {

	using namespace protocols;
	using namespace protocols::moves;

	// initialize core
	devel::init( argc, argv );

	MoverOP RosettaAntibody( new antibody_legacy::AntibodyModeler( ) );
	// protocols::jd2::JobDistributor::get_instance()->go( RosettaAntibody );
	protocols::jobdist::main_plain_mover( *RosettaAntibody );
	 } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

