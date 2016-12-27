// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

// libRosetta headers
#include <devel/init.hh>
#include <protocols/recces/sampler/MC_OneTorsion.hh>

// C++ headers
#include <iostream>
#include <fstream>

// Exception handling

using namespace core;
using namespace protocols::stepwise::sampler;

//////////////////////////////////////////////////////////////////////////////
void
run()
{
	core::id::TorsionID tor_id;
	MC_OneTorsion sampler(tor_id, 0);
	sampler.set_angle_range( -40, 40 );
	sampler.init();
	Size const n_step( 999999 );
	std::ofstream myfile;
	myfile.open( "angles.txt" );
	for ( Size i = 0; i != n_step; ++i ) {
		++sampler;
		sampler.update();
		myfile << sampler.active_angle() << std::endl;
	}
	myfile.close();
}
//////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	using namespace core;
	devel::init ( argc, argv );
	run();
  return 0;
}

