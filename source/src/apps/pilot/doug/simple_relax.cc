// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file demo/doug/rotamer_prediction_benchmark.cc
/// @brief this demo in conjuntion with a external script determins the percent of rotamers correctly predictied during a repack
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/Tracer.hh>

// utility headers
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>

// c++ headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <protocols/relax/ClassicRelax.hh>

#include <utility/excn/Exceptions.hh>


// namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;


void simple_relax( std::string pdb_filename, ScoreFunctionOP scorefxn );

int
main( int argc, char * argv [] )
{
    try {
	// initialize
	devel::init(argc, argv);

	// concatenate -s and -l flags together to get total list of PDB files
	// (This was taken from Ian's early job distributor, thanks Ian)
	std::vector< FileName > pdb_file_names;
	if ( option[ in::file::s ].active() )
		pdb_file_names = option[ in::file::s ]().vector(); // make a copy (-s)
	std::vector< FileName > list_file_names;
	if ( option[ in::file::l ].active() )
		list_file_names = option[ in::file::l ]().vector(); // make a copy (-l)

	for(std::vector< FileName >::iterator i = list_file_names.begin(), i_end = list_file_names.end(); i != i_end; ++i) {
		std::string filename( i->name() );
		std::ifstream data( filename.c_str() );
		if ( !data.good() ) {
			utility_exit_with_message( "Unable to open file: " + filename + '\n' );
		}
		std::string line;
		while( getline(data, line) ) {
			pdb_file_names.push_back( FileName(line) );
		}
		data.close();
	}

	// create score function
	ScoreFunctionOP scfxn( get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS ));

	// run simple relax for each name in list
	for(std::vector< FileName >::iterator i = pdb_file_names.begin(), i_end = pdb_file_names.end(); i != i_end; ++i) {
		simple_relax( i->name(), scfxn );
	}
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;
}

void simple_relax( std::string pdb_filename, ScoreFunctionOP scorefxn )
{
	// read in pose
	Pose pose;
	core::import_pose::pose_from_pdb( pose, pdb_filename );

	// create relax mover
	protocols::relax::ClassicRelax relax( scorefxn );

	// calculate orig score
	Energy orig_score = ( *scorefxn )( pose );

	// apply relax mover to pose
	std::cout << "Relaxing...";
	relax.apply(pose);
	std::cout << "Done!!!" << std::endl;

	// calculate score after relax
	Energy relax_score = ( *scorefxn )( pose );

	// print results
	std::cout << "Score before: " << orig_score << "\tScore after: " << relax_score << std::endl;

	// output pdb file
	io::pdb::dump_pdb( pose, pdb_filename + "_relax.pdb" );

	// ADD HBOND ANALYSIS CODE HERE
	// ADD HBOND ANALYSIS CODE HERE
	// ADD HBOND ANALYSIS CODE HERE
}
