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

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>


// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/util.hh>
#include <protocols/simple_moves/ShakeStructureMover.hh>
#include <utility/excn/Exceptions.hh>


using basic::Error;
using basic::Warning;


using namespace core;
using namespace core::fragment;
using namespace protocols;
using namespace utility;

using utility::vector1;

using core::io::pdb::dump_pdb;


namespace protocols {
namespace moves {


} // moves
} // protocols


///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		devel::init(argc, argv);
		using namespace protocols::moves;
		core::scoring::ScoreFunctionOP s( new core::scoring::ScoreFunction());
		protocols::moves::ShakeStructureMover *ssm = new ShakeStructureMover( );
		ssm->set_sc_min(true); //sc min after perturbing backbone
		ssm->set_nrounds( (500) );
		ssm->set_mc_temperature(  3 );
		MoverOP mover = ssm;
		protocols::jd2::JobDistributor::get_instance()->go( mover );
	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

