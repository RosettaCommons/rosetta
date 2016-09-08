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
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/util.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <core/pose/Pose.hh>


#include <basic/options/option.hh>
#include <protocols/viewer/viewers.hh>

#include <core/init/init.hh>

#include <utility/vector1.hh>

#include <ObjexxFCL/string.functions.hh>


// C++ headers
#include <iostream>
#include <string>


// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>

#include <utility/excn/Exceptions.hh>

using namespace core;
using namespace protocols;
using utility::vector1;


OPT_KEY( Boolean, remove_variant_cutpoint_atoms )

///////////////////////////////////////////////////////////////////////////////////////////
// Basically stolen from James' protein silent file extractor.
void
extract_pdbs_test()
{
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::io::silent;

	// setup residue types
	core::chemical::ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	// configure silent-file data object
	core::io::silent::SilentFileData silent_file_data;
	std::string infile  = option[ in::file::silent ][1];

	if ( option[ in::file::silent ].user() ) {
		silent_file_data.read_file( infile );
	} else {
		utility_exit_with_message( "Error: can't get any structures! Use -in::file::silent <silentfile>" );
	}

	core::pose::Pose pose;


	bool use_tags = false;
	std::set< std::string > desired_tags;
	if ( option[ in::file::tags ].active() ) {
		use_tags = true;
		desired_tags.insert( option[ in::file::tags ]().begin(), option[ in::file::tags ]().end() );
	}

	for ( core::io::silent::SilentFileData::iterator iter = silent_file_data.begin(), end = silent_file_data.end(); iter != end; ++iter ) {

		std::string const tag = iter->decoy_tag();

		if ( use_tags && ( !desired_tags.count( tag ) ) ) continue;

		std::cout << "Extracting: " << tag << std::endl;

		iter->fill_pose( pose, *rsd_set );

		//Fang: This outputing will fail if the pose contains virtual residue.
		//std::cout << "debug_rmsd(" << tag << ") = " << iter->get_debug_rmsd() << " over " << pose.size() << " residues... \n";


		if ( option[ remove_variant_cutpoint_atoms ]()==true ) {
			for ( Size n = 1; n <= pose.size(); n++  ) {
				pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_LOWER, n );
				pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_UPPER, n );
			}
		}

		std::string out_prefix = "";

		if ( option[out::prefix].user() ) out_prefix = option[out::prefix]();

		pose.dump_pdb( out_prefix + tag + ".pdb" );

	}

}
///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	extract_pdbs_test();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;

		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -in::file::silent <silent file> -tags <model tag like S_00001> " << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		NEW_OPT( remove_variant_cutpoint_atoms , "remove_variant_cutpoint_atoms", false );
		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		core::init::init(argc, argv);

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////

		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

