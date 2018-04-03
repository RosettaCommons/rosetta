// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
// //
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
// /// @file
// /// @brief  test application
// /// @author Melanie Aprahamian

#include <iostream>
#include <string>
#include <protocols/jd2/util.hh>
#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Main program
int main( int argc, char * argv [] )
{
	try
{

		// Initialize Rosetta
		devel::init( argc, argv );

		// Namespaces
		using namespace core;
		using namespace core::scoring;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace pose;

		// check if input PDB file is provided from command line with option -s
		if ( !option[in::file::s].user() ) {
			// exit if no PDB file is found
			utility_exit_with_message("Input PDB file not found");
		}

		//initialize pose
		pose::Pose p;

		//Get PDB file from command line option -s
		std::string pdb_file = option[in::file::s]()[1];

		// load pdb file into pose
		import_pose::pose_from_file(p,pdb_file);

		// initialize score function; score_fxn is a POINTER
		scoring::ScoreFunctionOP score_fxn = get_score_function();

		// score pose
		//(*score_fxn)(p);

		// output score
		score_fxn->show(p);

	}

catch ( utility::excn::EXCN_Base const & e )
{
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
	return 0;
}




