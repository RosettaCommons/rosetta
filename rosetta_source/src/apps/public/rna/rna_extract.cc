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


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED
#include <core/chemical/ChemicalManager.hh>

//#include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/rna/RNA_Util.hh>
// AUTO-REMOVED #include <core/scoring/rna/RNA_ScoringInfo.hh>

// AUTO-REMOVED #include <protocols/rna/RNA_ProtocolUtil.hh>

// AUTO-REMOVED #include <core/sequence/util.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
// AUTO-REMOVED #include <basic/options/option_macros.hh>

#include <protocols/viewer/viewers.hh>

#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <basic/basic.hh>
// AUTO-REMOVED #include <basic/database/open.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>

// AUTO-REMOVED #include <numeric/conversions.hh>
// AUTO-REMOVED #include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
//#include <cstdlib>
// AUTO-REMOVED #include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/io/mpistream.hh>



using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;

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
	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::RNA );

	// configure silent-file data object
	core::io::silent::SilentFileData silent_file_data;
	std::string infile  = option[ in::file::silent ][1];

	if ( option[ in::file::silent ].user() ) {
		silent_file_data.read_file( infile );
	} else {
		utility_exit_with_message( "Error: can't get any structures! Use -in::file::silent <silentfile>" );
	}

	core::pose::Pose pose;

	pose::Pose ideal_pose;
	//	bool const use_input_pose = option[ in::file::s ].active();
	//	if ( use_input_pose ) {
	//		std::string ideal_pdb_file  = option[ in::file::s ][1];
	//		core::import_pose::pose_from_pdb( ideal_pose, *rsd_set, ideal_pdb_file );
	//		protocols::rna::ensure_phosphate_nomenclature_matches_mini( ideal_pose );
	//	}

	bool use_tags = false;
	std::set< std::string > desired_tags;
	if( option[ in::file::tags ].active() ) {
		use_tags = true;
		desired_tags.insert( option[ in::file::tags ]().begin(), option[ in::file::tags ]().end() );
	}

	for ( core::io::silent::SilentFileData::iterator iter = silent_file_data.begin(), end = silent_file_data.end(); iter != end; ++iter ) {

		std::string const tag = iter->decoy_tag();

		if (use_tags && ( !desired_tags.count( tag ) ) ) continue;

		std::cout << "Extracting: " << tag << std::endl;

		//		if (use_input_pose) pose = ideal_pose;
		iter->fill_pose( pose, *rsd_set );

		std::cout << "debug_rmsd(" << tag << ") = " << iter->get_debug_rmsd() << " over " << pose.total_residue() << " residues... \n";
		pose.dump_pdb( tag + ".pdb" );

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
	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );

}
