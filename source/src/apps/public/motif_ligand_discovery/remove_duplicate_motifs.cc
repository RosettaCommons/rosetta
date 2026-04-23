// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief mdsmith, sthyme


#include <devel/init.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/io/pdb/pdb_writer.hh> // pose_from_pdb
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/dna/setup.hh>
#include <protocols/dna/RestrictDesignToProteinDNAInterface.hh>
#include <protocols/motifs/LigandMotifSearch.hh>
#include <protocols/motifs/LigandMotifSearch.fwd.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <protocols/motifs/Motif.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//#include <protocols/motifs/MotifLigandPacker.hh>
#include <protocols/dna/PDBOutput.hh>
#include <protocols/dna/util.hh>
//#include <protocols/dna/DnaInterfacePacker.hh>
//#include <protocols/dna/DnaInterfaceFinder.hh>
#include <protocols/motifs/motif_utils.hh>


#include <basic/prof.hh>
#include <basic/Tracer.hh>

#include <core/import_pose/import_pose.hh> // Need since refactor

// Utility Headers
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh> // file_exists
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>

#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>


// c++ headers
#include <fstream>
#include <iostream>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>

using namespace protocols::dna;
using utility::vector1;

using namespace core;
using namespace basic;
using namespace chemical;
using namespace pack;
using namespace task;
using namespace scoring;

////////////////////////////////////////////////////////////////////////////////
//This is a small app that serves to read in a motif library and remove any motifs that are functionally duplicates
//Functional duplicates are defined as motifs where all atom types match, and the distance and theta of the motifs are within a given cutoff
//non-duplicate motifs are returned in an outputted motifs library

//declare tracer
static basic::Tracer ms_tr( "apps.public.remove_duplicate_motifs" );

void
remove_duplicate_motifs(
	//core::pose::Pose & pose,
	//std::string & pdb_prefix,
	//protocols::motifs::MotifLibrary & motifs,
	//utility::io::ozstream & motif_output_file
)
{
	using utility::file::file_exists;

	DnaDesignDefOPs build_position_defs;

	using namespace core;
	using namespace ObjexxFCL;
	using namespace pose;
	using namespace chemical;
	using namespace scoring;
	using utility::vector1;

	//utility::vector1< utility::vector1< utility::vector1< core::Size > > > motif_indices_list; //Now, instead of having a vector of triplet vectors containing a size for the atom number, change size to vector of size to contain atom number and atomtype size
	core::Size ligand_marker = 1;
	utility::vector1< std::string > ligand_atom_names;

	protocols::motifs::MotifLibrary motifs( protocols::motifs::get_LigandMotifLibrary_user(ligand_marker, ligand_atom_names));
	protocols::motifs::MotifCOPs motif_library = motifs.library();

	//print initial size of motif library that has been read in
	ms_tr << "Size of read in motif library is " << motif_library.size() << std::endl;

	//Make new motif library that will have no duplicates
	protocols::motifs::MotifLibrary motif_library_no_duplicates;
	//core::Real duplicate_motif_cutoff = basic::options::option[ basic::options::OptionKeys::motifs::duplicate_motif_cutoff ](); // Default 0.2
	//iterate through motif library to search for duplicates

	//obtain thresholds for determining when a motif is considered a duplicate for angle and distance between motifs
	Real dist_threshold( basic::options::option[ basic::options::OptionKeys::motifs::duplicate_dist_cutoff ]  );
	Real angl_threshold( basic::options::option[ basic::options::OptionKeys::motifs::duplicate_angle_cutoff ] );

	for ( auto single_motif : motif_library ) {

		//single_motif->place_residue( *protres, *res );
		bool broke( false );

		for ( auto motifcop : motif_library_no_duplicates ) {

			//check if atoms match and ensure same residue names on both sides
			if ( motifcop->motif_atom_match_strict(*single_motif) == false ) continue;

			//motifcop->place_residue( *protres, *res2 );


			core::Real distance = 0;
			core::Real theta = 0;



			jump_distance(single_motif->forward_jump(), motifcop->forward_jump(), distance, theta);
			//ms_tr << "Comparing to motif " << motifcop->remark() << std::endl;
			//ms_tr << "Distance: " << distance << std::endl;
			//ms_tr << "Theta: " << theta << std::endl;

			//core::Real rmsdtest = 0.2;//core::scoring::automorphic_rmsd( *res, *res2, false );
			if ( distance < dist_threshold && theta < angl_threshold ) {
				ms_tr << "Skipping motif with distance difference: " << distance << "  and angle difference: "  << theta << std::endl;
				if ( single_motif->has_remark() ) {
					ms_tr << "Skipped motif remark is: " << single_motif->remark() << std::endl;
				}
				if ( motifcop->has_remark() ) {
					ms_tr << "Skipped motif matches motif with remark: " << motifcop->remark() << std::endl;
				}

				broke = true;
				break;
			}
		}
		if ( !broke ) {
			ms_tr << "Adding motif " << single_motif->remark() << std::endl;
			//final_contacts.push_back( contacts_strength_filter[ic] );
			motif_library_no_duplicates.add_to_library( *single_motif );
		}
	}

	//run through built motif_library_no_duplicates and write to a new motifs file
	ms_tr << "Adding motifs from a previously made library" << std::endl;
	std::string MotifLibraryFileName = "MotifLibrary.motifs";
	if ( basic::options::option[ basic::options::OptionKeys::motifs::output_file ].user() ) {
		MotifLibraryFileName = basic::options::option[ basic::options::OptionKeys::motifs::output_file ]();
	}
	utility::io::ozstream motif_output_file( MotifLibraryFileName );
	int final_motif_counter = 0;
	for ( auto motifcop : motif_library_no_duplicates ) {
		//motifs.add_to_library( *motifcop );
		motif_output_file << *motifcop;
		++final_motif_counter;
	}
	ms_tr << "Final number of motifs when filtering for duplicates is: " << final_motif_counter << std::endl;
}



int
main( int argc, char * argv [] )
{
	try {
		using namespace options;
		using namespace OptionKeys;

		devel::init( argc, argv ); // reading options--name should be more descriptive

		basic::prof_reset();

		ms_tr << "Getting input motif file." << std::endl;

		if ( option[ OptionKeys::motifs::list_motifs ].user() || option[ OptionKeys::motifs::motif_filename ].user() ) {
			ms_tr <<  "Removing duplicates." << std::endl;
			remove_duplicate_motifs();
			ms_tr << "SUCCESSFUL COMPLETION" << std::endl;
		} else {
			ms_tr << "Neither the list_motifs nor motif_filename flags were used, and so no motif list was provided." << std::endl;
		}
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
