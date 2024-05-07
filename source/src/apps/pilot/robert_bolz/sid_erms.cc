// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <iostream>
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

//local options
basic::options::StringOptionKey const expERMS( "expERMS" );
basic::options::StringOptionKey const predERMS( "predERMS" );

int main( int argc, char ** argv ) {
	basic::options::option.add( expERMS, "Experimental ERMS from SID" ).def("");
	basic::options::option.add( predERMS, "predicted ERMS from SID" ).def("");
	devel::init( argc, argv );
	utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
	if ( filenames.size() > 0 ) {
		std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;	
	} else {
		std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
		return 1;
	}

	core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );

	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
	
	core::Real score = sfxn->score( *mypose );
	std::cout << "Score is: " << score << std::endl;


	std::string exp_ERMS = basic::options::option[ expERMS ].value();
	std::string pred_ERMS = basic::options::option[ predERMS ].value();
	//std::istringstream ExpERMS; //= utility::io::izstream::open(exp_ERMS)
	//std::istringstream PredERMS; //= utility::io::izstream::open(pred_ERMS)
	std::string ExpERMS;
	std::string PredERMS;
	std::fstream(exp_ERMS, std::ios::in ) >> ExpERMS;
	std::fstream(pred_ERMS, std::ios::in ) >> PredERMS;
	std::cout << "exp ERMS is: " << ExpERMS << std::endl;
	std::cout << "pred ERMS is: " << PredERMS << std::endl;

	return 0;
}