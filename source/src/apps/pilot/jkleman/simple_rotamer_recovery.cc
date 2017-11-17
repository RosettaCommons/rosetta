// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief   Computes rotamer recovery between an input file and a native
/// @details last Modified: 08/26/15
/// @author  JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Package Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>
#include <core/pack/util.hh>

// utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <cstdlib>
#include <fstream>

using basic::Error;
using basic::Warning;
using std::ofstream;

static basic::Tracer TR( "apps.pilot.jkleman.simple_rotamer_recovery" );

////////////////////////////////////////////////////////////////////////////////

void simple_rotamer_recovery(){

	using namespace core;
	using namespace core::pack;
	using namespace core::pose;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// read in PDB or list of PDBs
	utility::vector1< std::string > filenames;
	utility::vector1< Pose > poses;
	Pose pose, native;
	std::string native_file;

	// if single PDB
	if ( option[OptionKeys::in::file::s].user() ) {
		std::string filename = option[ OptionKeys::in::file::s ].value_string();
		core::import_pose::pose_from_file( pose, filename , core::import_pose::PDB_file);
		poses.push_back( pose );
		filenames.push_back( filename );
	} else if ( option[OptionKeys::in::file::l].user() ) {
		// if list of PDBs

		// read filenames
		filenames = option[OptionKeys::in::file::l]();

		// add pose to vector of poses
		for ( core::Size i = 1; i <= filenames.size(); ++i ) {
			Pose pdb;
			core::import_pose::pose_from_file( pdb, filenames[ i ] , core::import_pose::PDB_file);
			poses.push_back( pdb );
		}
	} else {
		utility_exit_with_message( "Please provide a PDB file with the option -in:file:s or -in:file:l. Quitting." );
	}

	// read in native
	if ( option[ OptionKeys::in::file::native ].user() ) {
		native_file = option[ OptionKeys::in::file::native ].value_string();
		core::import_pose::pose_from_file( native, native_file , core::import_pose::PDB_file);
	} else {
		utility_exit_with_message( "Please provide a native PDB file with the option -in:file:native. Quitting." );
	}

	// open output file
	ofstream fout;
	fout.open( "rotamer_recovery.out", std::ios::out );
	if ( !fout.is_open() ) {
		TR << "Unable to open output file 'rotamer_recovery.out'." << std::endl;
		utility_exit();
	}

	// iterate through PDBs
	for ( core::Size i = 1; i <= poses.size(); ++i ) {

		// angle difference
		core::Real difference = 10.0;

		// compute rotamer recovery and residue rotamer recovery
		core::Real rot_rec( core::pack::rotamer_recovery( poses[ i ], native, difference ) );
		core::Real res_rot_rec( core::pack::residue_rotamer_recovery( poses[ i ], native, difference ) );

		// write output
		fout << "rotamer recovery between " << native_file << " and " << filenames[ i ] << " " << rot_rec << std::endl;
		fout << "residue rotamer recovery between " << native_file << " and " << filenames[ i ] << " " << res_rot_rec << std::endl;

	}
	fout.close();

} // simple rotamer recovery

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// MAIN ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {

		// initialize option system, RNG, and all factory-registrators
		devel::init(argc, argv);

		// call my function
		simple_rotamer_recovery();

	}
catch (utility::excn::Exception const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
}
