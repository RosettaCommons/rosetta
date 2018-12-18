// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file get_binding_site.cc
/// @author Noah Ollikainen
/// @author Roland A. Pache, PhD

// Core headers
#include <devel/init.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <core/pack/task/residue_selector/ClashBasedShellSelector.hh>
#include <basic/Tracer.hh>

// Option headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <iostream>
#include <fstream>


int
main( int argc, char * argv [] )
{
	try{

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::pack::task;

		// initialize Rosetta
		devel::init(argc, argv);
		core::pack::task::TaskFactoryOP main_task_factory( new TaskFactory );

		// read resfile and other parameters from command line
		main_task_factory->push_back( operation::TaskOperationCOP( new operation::InitializeFromCommandline ) );
		if ( option[ packing::resfile ].user() ) {
			main_task_factory->push_back( operation::TaskOperationCOP( new operation::ReadResfile ) );
		} else {
			std::cout << "WARNING: no resfile given in command line, restricting all residues to repacking" << std::endl;
			main_task_factory->push_back( operation::TaskOperationCOP( new operation::RestrictToRepacking() ) );
		}

		// setup pose and score function
		core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();
		core::pose::Pose starting_pose = core::pose::Pose();
		core::import_pose::pose_from_file( starting_pose, option[ in::file::s ]()[1] , core::import_pose::PDB_file);
		core::pack::task::PackerTaskOP task( main_task_factory->create_task_and_apply_taskoperations( starting_pose ) );

		// make a copy of the original pose
		core::pose::PoseOP pose_copy( new core::pose::Pose(starting_pose) );

		// using ClashBasedShellSelector to define repack shell
		// old implementation of ClashBasedShellSelector
		// core::pack::task::residue_selector::ClashBasedShellSelectorOP cbrss( new core::pack::task::residue_selector::ClashBasedShellSelector(task, score_fxn) );
		// new implementation of ClashBasedShellSelector
		core::pack::task::residue_selector::ClashBasedShellSelectorOP cbrss( new core::pack::task::residue_selector::ClashBasedShellSelector(task) );
		cbrss->set_include_focus(false);
		utility::vector1< bool > to_repack = cbrss->apply( *pose_copy );

		// print list of positions to repack
		std::cout << std::endl << "Clash-based repack shell:" << std::endl;
		core::Size num_positions_to_repack=0;
		for ( core::Size i = 1; i <= to_repack.size(); ++i ) {
			if ( to_repack[i] ) {
				std::cout << pose_copy->pdb_info()->number(i) << " "
					<< pose_copy->pdb_info()->chain(i) << " NATAA" << std::endl;
				++num_positions_to_repack;
			}
		}
		std::cout << std::endl << "For highlighting in pymol:" << std::endl;
		std::stringstream for_pymol;
		for ( core::Size i = 1; i <= to_repack.size(); ++i ) {
			if ( to_repack[i] ) {
				for_pymol << "chain " << pose_copy->pdb_info()->chain(i) << " and resi " << pose_copy->pdb_info()->number(i) << " or ";
			}
		}
		std::string for_pymol_string=for_pymol.str();
		std::cout << for_pymol_string.substr(0,for_pymol_string.size()-4) << std::endl;
		std::cout << std::endl << "Found " << num_positions_to_repack << " positions to repack" << std::endl;

		// append list of repack positions to resfile
		if ( option[ packing::resfile ].user() ) {
			std::string resfile_name=option[ packing::resfile ][1];
			std::cout << "Appending list of repack positions to " << resfile_name << std::endl;
			std::ofstream resfile( resfile_name.c_str(), std::ios::app );
			if ( resfile.is_open() ) {
				// store list of positions to repack
				for ( core::Size i = 1; i <= to_repack.size(); ++i ) {
					if ( to_repack[i] ) {
						resfile << std::endl << pose_copy->pdb_info()->number(i) << " "
							<< pose_copy->pdb_info()->chain(i) << " NATAA";
					}
				}
				resfile.close();
			} else std::cout << "ERROR: Unable to open resfile" << std::endl;
		}

		return 0;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
