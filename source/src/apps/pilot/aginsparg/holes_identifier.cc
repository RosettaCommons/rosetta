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
//#include <core/scoring/DockingScoreFunction.hh>
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

#include <protocols/simple_filters/HolesFilter.hh>
#include <protocols/protein_interface_design/filters/InterfaceHolesFilter.hh>

//file path navigation
//#include <experimental/filesystem>

//for data type debugging
#include <typeinfo>

using namespace protocols::dna;
using utility::vector1;
using utility::string_split;
using namespace core;
using namespace basic;
using namespace chemical;
using namespace pack;
using namespace task;
using namespace scoring;
//using namespace apps;

int
main( int argc, char * argv [] )
{
	try {
		//example run:
		// /data/project/thymelab/running_Rosetta/ari_work/Rosetta_Code_copy/main/source/bin/holes_identifier.linuxgccrelease -s /data/project/thymelab/running_Rosetta/ari_work/REU_shared_space/s_lib/scripts/ari_discovery_results_HIS253/score_-25/0/good_pdb.pdb -holes:dalphaball /data/project/thymelab/running_Rosetta/ari_work/Rosetta_Code_copy/main/source/external/DAlpahBall/DAlphaBall.gcc -extra_res_fa good_params.params

		using namespace options;
		using namespace OptionKeys;

		static basic::Tracer ms_tr( "Code tracing", basic::t_info );

		devel::init( argc, argv ); // reading options--name should be more descriptive

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


		//loop through inputted pdbs
		for ( Filenames::const_iterator filename( pdbnames.begin() );
				filename != pdbnames.end(); ++filename ) {
			if ( !utility::file::file_exists( *filename ) ) {
				continue;
			}

			pose::Pose my_pose;

			core::import_pose::pose_from_file( my_pose, *filename , core::import_pose::PDB_file);

			//make copy of pose in case apply modifies original
			pose::Pose copy_pose = my_pose;


			//make HolesFilter object
			protocols::simple_filters::HolesFilter holesfilter;

			holesfilter.set_exclude_bb_atoms(false);

			//run apply and compute on the pose and output it
			bool apply_result = holesfilter.apply(my_pose);

			core::Real compute_result = holesfilter.compute(copy_pose);

			ms_tr << *filename << ", " << apply_result << ", " << compute_result << std::endl;

			/*
			holesfilter.set_exclude_bb_atoms(true);

			bool apply_result2 = holesfilter.apply(my_pose);

			core::Real compute_result2 = holesfilter.compute(copy_pose);

			ms_tr << *filename << " with exclude bb atoms, " << apply_result2 << ", " << compute_result2 << std::endl;

			protocols::protein_interface_design::filters::InterfaceHolesFilter ihf;

			bool apply_result3 = ihf.apply(my_pose);

			core::Real compute_result3 = ihf.compute(copy_pose);

			ms_tr << *filename << " Inferface Holes Filter, " << apply_result3 << ", " << compute_result3 << std::endl;
			*/

			/*
			//run apply and compute on the pose and output it
			bool apply_result = protocols::simple_filters::HolesFilter::apply(my_pose);

			core::Real compute_result = protocols::simple_filters::HolesFilter::compute(copy_pose);

			ms_tr << filename << ", " << apply_result << ", " << compute_result << std::endl;
			*/

		}



	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
