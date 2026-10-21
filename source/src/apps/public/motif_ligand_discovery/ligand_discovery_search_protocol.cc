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

#include <core/import_pose/import_pose.hh> // Need since refactor

#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>

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
#include <string>
#include <queue>
#include <functional>
#include <tuple>

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
using utility::string_split;
using utility::vector1;

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
		bool check_for_bad_motifs = true;
		protocols::motifs::MotifLibrary motifs( protocols::motifs::get_LigandMotifLibrary_user(check_for_bad_motifs, ligand_atom_names) );
		protocols::motifs::MotifCOPs motifcops = motifs.library();
		ms_tr << "Collected motifs library with " << motifcops.size() << " motifs. " << std::endl;

		bool const use_input_ligand = option[ OptionKeys::motifs::input_ligand_residue ].user();
		std::string const input_ligand_start_step = option[ OptionKeys::motifs::input_ligand_start_step ];
		if ( !use_input_ligand && option[ OptionKeys::motifs::input_ligand_start_step ].user() ) {
			ms_tr.Error << "input_ligand_start_step requires input_ligand_residue." << std::endl;
			return -1;
		}
		if ( !use_input_ligand && option[ OptionKeys::motifs::input_ligand_skip_highresdock ].user() ) {
			ms_tr.Error << "input_ligand_skip_highresdock requires input_ligand_residue." << std::endl;
			return -1;
		}
		int input_ligand_pdb_number = 0;
		char input_ligand_chain = ' ';
		char input_ligand_insertion_code = ' ';
		if ( use_input_ligand ) {
			auto const residue_selector = option[ OptionKeys::motifs::input_ligand_residue ].resnum_and_chain();
			utility::vector1<int> const & selected_residue_numbers = std::get<0>(residue_selector);
			utility::vector1<char> const & selected_chains = std::get<1>(residue_selector);
			if ( selected_residue_numbers.size() != 1 || selected_chains.size() != 1 ) {
				ms_tr.Error << "input_ligand_residue must identify exactly one PDB residue, e.g. A:479." << std::endl;
				return -1;
			}
			input_ligand_pdb_number = selected_residue_numbers[1];
			input_ligand_chain = selected_chains[1];
			if ( option[ OptionKeys::motifs::input_ligand_insertion_code ].user() ) {
				std::string const insertion_code = option[ OptionKeys::motifs::input_ligand_insertion_code ];
				if ( insertion_code.size() != 1 ) {
					ms_tr.Error << "input_ligand_insertion_code must be exactly one character." << std::endl;
					return -1;
				}
				input_ligand_insertion_code = insertion_code[0];
			}
		} else if ( option[ OptionKeys::motifs::input_ligand_insertion_code ].user() ) {
			ms_tr.Error << "input_ligand_insertion_code requires input_ligand_residue." << std::endl;
			return -1;
		}

		//load in ligands

		std::string tag = "fa_standard";
		std::string directory = "";
		core::chemical::ResidueTypeCOPs all_rsd_types;
		core::chemical::GlobalResidueTypeSetOP input_residue_set;

		if ( option[ OptionKeys::motifs::params_directory_path ].user() ) {
			directory = option[ OptionKeys::motifs::params_directory_path ];
		}

		//uses params_directory_path flag
		if ( option[ OptionKeys::motifs::params_directory_path ].user() ) {
			input_residue_set = core::chemical::GlobalResidueTypeSetOP(new core::chemical::GlobalResidueTypeSet(tag, directory));
			core::chemical::ResidueTypeCOPs const & temp_all_rsd_types( core::chemical::ResidueTypeFinder( *input_residue_set ).get_all_possible_residue_types() );
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
			if ( use_input_ligand && input_residue_set != nullptr ) {
				core::import_pose::pose_from_file(pose_pointer_helper, *input_residue_set, *filename, false, core::import_pose::PDB_file);
			} else {
				core::import_pose::pose_from_file( pose_pointer_helper, *filename , core::import_pose::PDB_file);
			}
			pose::PoseOP pose(new pose::Pose(pose_pointer_helper));

			//get discovery position(s), load into Size vector
			//vector from options is int, but assignment should convert to size
			utility::vector1<core::Size> discovery_position = option[ OptionKeys::motifs::protein_discovery_locus ]();

			//declare LigandDiscoverySearch object
			ms_tr << "Making LigandDiscoverySearch object" << std::endl;
			//protocols::motifs::LigandDiscoverySearch lds(pose, motifcops, all_residues, discovery_position);
			LigandDiscoverySearch lds(pose, motifcops, all_residues, discovery_position);

			//optionally seed the lds with secondary residue values if the motifs::secondary_protein_discovery_locus flag was used
			if ( option[ OptionKeys::motifs::secondary_protein_discovery_locus ].user() ) {
				//set inputs as Size
				utility::vector1<core::Size> secondary_discovery_position = option[ OptionKeys::motifs::secondary_protein_discovery_locus ]();
				//use vector overload of setter
				lds.set_secondary_working_positions(secondary_discovery_position);
			}

			//debug option to easily spit out the empty receptor that has residue indices converted to Rosetta indexing 
			if ( option[ OptionKeys::motifs::dump_input_pdb_with_rosetta_residue_indices ].user() ) {
				core::io::pdb::dump_pdb(*pose, "empty_rosettaindices_" + pdbprefix + ".pdb");
			}

			//run discovery
			ms_tr << "Running discovery for " << pdbprefix << std::endl;
			if ( use_input_ligand ) {
				if ( pose->pdb_info() == nullptr ) {
					ms_tr.Error << "Input pose has no PDB numbering information; cannot locate the selected ligand." << std::endl;
					return -1;
				}
				core::Size ligand_residue_position = 0;
				core::Size matching_residue_count = 0;
				for ( core::Size residue_position = 1; residue_position <= pose->size(); ++residue_position ) {
					if ( pose->pdb_info()->chain(residue_position) == input_ligand_chain &&
							pose->pdb_info()->number(residue_position) == input_ligand_pdb_number &&
							pose->pdb_info()->icode(residue_position) == input_ligand_insertion_code ) {
						ligand_residue_position = residue_position;
						++matching_residue_count;
					}
				}
				if ( matching_residue_count != 1 ) {
					ms_tr.Error << "Expected exactly one PDB residue matching chain " << input_ligand_chain
						<< ", residue " << input_ligand_pdb_number << ", insertion code '"
						<< input_ligand_insertion_code << "'; found " << matching_residue_count << "." << std::endl;
					return -1;
				}
				if ( !pose->residue(ligand_residue_position).is_ligand() ) {
					ms_tr.Error << "The selected PDB residue is not recognized as a ligand: "
						<< pose->residue(ligand_residue_position).name() << std::endl;
					return -1;
				}
				lds.discover(pdbprefix, ligand_residue_position, input_ligand_start_step);
			} else {
				lds.discover(pdbprefix);
			}
		}
	}
catch (utility::excn::Exception const & e ) {
	e.display();
	return -1;
}
	return 0;
}
//main/source/src/protocols_a.5.src.settings
