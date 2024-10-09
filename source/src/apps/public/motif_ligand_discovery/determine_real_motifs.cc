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
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <sstream>
#include <string>

#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>

#include <protocols/motifs/IdentifyLigandMotifs.fwd.hh>
#include <protocols/motifs/IdentifyLigandMotifs.hh>
#include <protocols/motifs/Motif.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <protocols/motifs/motif_utils.hh>

#include <core/io/pdb/pdb_writer.hh> // pose_from_pdb
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

#include <core/import_pose/import_pose.hh> // Need since refactor

#include <basic/options/option.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh> // file_exists
#include <utility/vector1.hh>

// Time profiling header
#include <time.h>

using namespace core;
//using namespace ObjexxFCL;
using namespace pose;
//using namespace chemical;
using namespace scoring;
//using namespace optimization;
using namespace basic;
using namespace options;
using namespace OptionKeys;

using utility::string_split;
using utility::vector1;

int
main( int argc, char * argv [] )
{
	try {
		static basic::Tracer ms_tr( "determine_real_motifs", basic::t_info );



		devel::init( argc, argv ); // reading options--name should be more descriptive

		//load in pdb(s)

		ms_tr << "Making Pose of input pdb file(s)" << std::endl;

		typedef vector1< utility::file::FileName > Filenames;
		Filenames pdbnames;

		//borrowing code from motif_ligand_packer_design.cc to take in a single or multiple pdb files (using -l or -s flags), and be able to iterate discovery upon all files
		if ( option[ in::file::l ].user() ) {
			Filenames listnames( option[ in::file::l ]().vector() );
			for ( Filenames::const_iterator filename( listnames.begin() );
					filename != listnames.end(); ++filename ) {
				std::ifstream list( (*filename).name().c_str() );
				while ( list ) {
					std::string pdbname;
					list >> pdbname;
					pdbnames.push_back( pdbname );
				}
			}

		} else if ( option[ in::file::s ].user() ) {
			pdbnames = option[ in::file::s ]().vector();

		} else {
			std::cerr << "No files given: Use either -file:s or -file:l "
				<< "to designate a single pdb or a list of pdbs"
				<< std::endl;
		}

		//read in the compare library of motifs (presummably real motifs)
		//create a vector that only stores the identity of the atoms  in the ligand (not their index  value); this vector will be empty, but is necessary to call wanted versions of overloaded  function to loading in the library
		//version if we don't use this empty vector has potential to terminate the read of motifs early if a single bad motif is  read in, this version ditches the bad motif and keeps reading
		utility::vector1< std::string > ligand_atom_names;
		bool check_for_bad_motifs = true;
		protocols::motifs::MotifLibrary real_motifs( protocols::motifs::get_LigandMotifLibrary_user(check_for_bad_motifs, ligand_atom_names) );
		protocols::motifs::MotifCOPs real_motifcops = real_motifs.library();
		ms_tr << "Collected motifs library with " << real_motifcops.size() << " motifs. " << std::endl;

		//hash the library by residue and atoms
		std::map<protocols::motifs::motif_atoms,protocols::motifs::MotifCOPs> hashed_real_library;
		protocols::motifs::hash_motif_library_into_map(real_motifs,hashed_real_library);
		
		//create and use IdentifyLigandMotifs
		IdentifyLigandMotifs ilm;

		//iterate through all pdb files and run the real motif check
		for ( Filenames::const_iterator filename( pdbnames.begin() ); filename != pdbnames.end(); ++filename ) {

			if ( !utility::file::file_exists( *filename ) ) {
				continue;
			}

			//get pdb name prefix for outputs:
			//strip path and file type
			std::string pdbprefix( string_split( string_split( *filename, '/' ).back(), '.' ).front() );

			//import pose from pdb
			pose::Pose pose_pointer_helper;
			core::import_pose::pose_from_file( pose_pointer_helper, *filename , core::import_pose::PDB_file);
			pose::PoseOP pose(new pose::Pose(pose_pointer_helper));

			//run evaluate_motifs_of_pose() on the pose
			ilm.evaluate_motifs_of_pose(pose, hashed_real_library, pdbprefix);

			//dump the pose to a pdb file
			core::io::pdb::dump_pdb(*pose, pdbprefix + ".pdb");
		}

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
