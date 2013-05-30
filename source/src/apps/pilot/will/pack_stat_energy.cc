// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


#include <core/types.hh>
#include <devel/init.hh>

#include <basic/options/option.hh>

#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <utility/io/izstream.hh>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>






void test( std::string fname )
{
	using namespace core;
	using namespace pose;
	using namespace scoring;

	Pose pose;
	core::import_pose::pose_from_pdb(pose,fname);

  ScoreFunction sf( *(getScoreFunctionLegacy( PRE_TALARIS_2013_STANDARD_WTS )) );

  ScoreFunction sfpack;
	sfpack.set_weight(pack_stat,1.0);

	std::cerr << fname << " " << sf(pose) << std::endl;
	std::cerr << fname << " " << sfpack(pose) << std::endl;

}

int
main (int argc, char *argv[])
{

	try {



	devel::init( argc, argv );

  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace utility;

  // test_io();

	// test_sasa_dots();

	if( option[ in::file::s ].user() ) {
  	vector1<file::FileName> files( option[ in::file::s ]() );
  	for( size_t i = 1; i <= files.size(); ++i ) {
    	test( files[i] );
  	}
	} else if( option[ in::file::l ].user() ) {
  	vector1<file::FileName> files( option[ in::file::l ]() );
  	for( size_t i = 1; i <= files.size(); ++i ) {
			utility::io::izstream list( files[i] );
			std::string fname;
			while( list >> fname ) {
				// std::cerr << "'" << fname << "'" << std::endl;
    		test( fname );
			}
  	}
	}
	return 0;



	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}
