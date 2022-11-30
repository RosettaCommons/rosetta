// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief aginsparg, ipatel, sthyme


#include <devel/init.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/io/pdb/pdb_writer.hh> // pose_from_pdb
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/dna/setup.hh>
#include <protocols/dna/RestrictDesignToProteinDNAInterface.hh>
#include <protocols/motifs/LigandMotifSearch.hh>
#include <protocols/motifs/LigandMotifSearch.fwd.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ChemicalManager.hh>
#include <protocols/motifs/Motif.hh>
#include <core/chemical/residue_io.hh>
#include <protocols/motifs/MotifHit.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <protocols/motifs/BuildPosition.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>

#include <core/pose/PDBInfo.hh>

#include <core/chemical/GlobalResidueTypeSet.hh>

//#include <protocols/motifs/MotifLigandPacker.hh>
#include <protocols/dna/PDBOutput.hh>
#include <protocols/dna/util.hh>
//#include <protocols/dna/DnaInterfacePacker.hh>
//#include <protocols/dna/DnaInterfaceFinder.hh>
#include <protocols/motifs/motif_utils.hh>
using namespace protocols::dna;

#include <basic/prof.hh>
#include <basic/Tracer.hh>
static basic::Tracer TR( "apps.pilot.motif_dna_packer_design" );

#include <core/import_pose/import_pose.hh> // Need since refactor

#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>

// Utility Headers
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh> // file_exists
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
using utility::vector1;
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
using utility::string_split;

// c++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <queue>
#include <functional>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>

//this may solve atom clash issue
#include <core/pose/xyzStripeHashPose.hh>
#include <numeric/geometry/hashing/xyzStripeHash.hh>
#include <numeric/geometry/hashing/xyzStripeHash.fwd.hh>

#include <numeric/xyzVector.hh>
#include <protocols/ligand_docking/InterfaceScoreCalculator.hh>
#include <protocols/ligand_docking/ligand_scores.hh>

#include <protocols/qsar/scoring_grid/AtrGrid.hh>
#include <protocols/qsar/scoring_grid/RepGrid.hh>
#include <protocols/qsar/scoring_grid/VdwGrid.hh>
#include <protocols/qsar/scoring_grid/HbaGrid.hh>
#include <protocols/qsar/scoring_grid/HbdGrid.hh>
#include <protocols/qsar/scoring_grid/GridSet.hh>

#include <protocols/ligand_docking/MoveMapBuilder.hh>
#include <protocols/ligand_docking/LigandArea.hh>
#include <protocols/ligand_docking/InterfaceBuilder.hh>
#include <protocols/ligand_docking/MoveMapBuilder.fwd.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/ligand_docking/LigandBaseProtocol.hh>
#include <protocols/ligand_docking/HighResDocker.hh>
#include <protocols/ligand_docking/LigandDockProtocol.hh>
#include <protocols/ligand_docking/ligand_dock_impl.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/chemical/PoseResidueTypeSet.hh>

#include <ObjexxFCL/FArray1D.hh>

// Time profiling header
#include <time.h>

#include <protocols/motifs/LigandDiscoverySearch.hh>

//for data type debugging
#include <typeinfo>

using namespace core;
using namespace basic;
using namespace chemical;
using namespace pack;
using namespace task;
using namespace scoring;

////////////////////////////////////////////////////////////////////////////////



int
main( int argc, char * argv [] )
{
	try
{

		static basic::Tracer ms_tr( "LDS_app", basic::t_info );

		using namespace options;
		using namespace OptionKeys;

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

		//load in whole motif library
		//create a vector that only stores the identity of the atoms  in the ligand (not their index  value); this vector will be empty, but is necessary to call wanted versions of overloaded  function to loading in the library
		//version if we don't use this empty vector has potential to terminate the read of motifs early if a single bad motif is  read in, this version ditches the bad motif and keeps reading
		utility::vector1< std::string > ligand_atom_names;
		core::Size const ligand_marker = 0;
		protocols::motifs::MotifLibrary motifs( protocols::motifs::get_LigandMotifLibrary_user(ligand_marker, ligand_atom_names) );
		protocols::motifs::MotifCOPs motifcops = motifs.library();

		//load in ligands

		std::string tag = "fa_standard";
		std::string directory = "";
		core::chemical::ResidueTypeCOPs all_rsd_types;

		if ( option[ OptionKeys::motifs::params_directory_path ].user() ) {
			directory = option[ OptionKeys::motifs::params_directory_path ];
		}

		//uses params_directory_path flag
		if ( option[ OptionKeys::motifs::params_directory_path ].user() ) {
			core::chemical::GlobalResidueTypeSet residue_set(tag, directory);
			core::chemical::ResidueTypeCOPs const & temp_all_rsd_types( core::chemical::ResidueTypeFinder( residue_set ).get_all_possible_residue_types() );
			all_rsd_types = temp_all_rsd_types;
		}

		//convert all_rsd_types into a vector of residue (from residue.cc) data type
		utility::vector1<core::conformation::ResidueOP> all_residues;

		//run through all_rsd_types and convert to residue
		for ( auto ref : all_rsd_types ) {
			core::chemical::ResidueTypeCOP ligres(ref);
			//convert ligres to be a ResidueOP type
			core::conformation::ResidueOP ligresOP(new core::conformation::Residue (ligres, true));

			all_residues.push_back(ligresOP);
		}

		//iterate through all pdb files to run discovery on
		for ( Filenames::const_iterator filename( pdbnames.begin() );
				filename != pdbnames.end(); ++filename ) {

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

			//get discovery position
			core::Size discovery_position = option[ OptionKeys::motifs::protein_discovery_locus ];

			//declare LigandDiscoverySearch object
			ms_tr << "Making LigandDiscoverySearch object" << std::endl;
			//protocols::motifs::LigandDiscoverySearch lds(pose, motifcops, all_residues, discovery_position);
			LigandDiscoverySearch lds(pose, motifcops, all_residues, discovery_position);
			//run discovery
			ms_tr << "Running discovery for " << pdbprefix << std::endl;
			lds.discover(pdbprefix);
		}
	}
catch (utility::excn::Exception const & e ) {
	e.display();
	return -1;
}
	return 0;
}
//main/source/src/protocols_a.5.src.settings
