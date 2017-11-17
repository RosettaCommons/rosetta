// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
//
/// @file   FlexPepDock.cc
//
/// @brief Main file for applying flexible peptide docking protocol
/// @author Barak Raveh
/// @date August 05, 2008

//#define GL_GRAPHICS

// Mini-Rosetta headers
#include <devel/FlexPepDocking/FlexPepDockingProtocol.hh>

#include <numeric/random/random.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/util.hh>//option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/basic.hh>
#include <basic/Tracer.hh>

//#include <protocols/jobdist/standard_mains.hh>
//#include <protocols/moves/Mover.hh>

#include <ObjexxFCL/format.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>

using basic::Error;
using basic::Warning;


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		using namespace core;

		devel::init(argc, argv);

		pose::Pose pose;
		core::import_pose::pose_from_file( pose, basic::options::start_file() , core::import_pose::PDB_file);
		using namespace core::scoring;
		// standard packer wts
		std::cout << "Standard packer weights without fa_rep term:" << std::endl;
		ScoreFunctionOP scorefxn( get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS ) );
		scorefxn->set_weight( fa_rep, 0.0 );
		(*scorefxn)(pose);
		scorefxn->show(std::cout,pose);

		//  std::cout << *scorefxn;

		/// soft rep packer wts
		std::cout << "Soft repulsive weights" << std::endl;
		ScoreFunctionOP scorefxn2( ScoreFunctionFactory::create_score_function( SOFT_REP_WTS ) );
		(*scorefxn2)(pose);

		std::cout << *scorefxn2;
		scorefxn2->show(std::cout,pose);

		/// score12 w/ std packer wts
		std::cout << "Score12 with standard packer weights:" << std::endl;
		ScoreFunctionOP score12( get_score_function() );
		(*score12)(pose);

		std::cout << *score12;
		score12->show(std::cout,pose);


	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
