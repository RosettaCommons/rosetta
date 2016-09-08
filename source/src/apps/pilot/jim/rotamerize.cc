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

#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/util.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <protocols/simple_moves/RotamerizeMover.hh>

#include <devel/init.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <set>
#include <vector>
#include <list>
#include <iterator>

using namespace std;
using namespace core;
using namespace core::id;
using namespace core::scoring;
using namespace protocols::simple_moves;
using namespace utility;

void
rotamerize( std::string filename, pose::Pose & pose )
{
	pose::Pose start_pose;
	start_pose = pose;

	std::string weights( "talaris2013" );

	ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( weights ) );;

	Energy start_score = (*score_fxn)( pose );
//	pose.energies().show( std::cout );
//	std::cout << "Starting energy is " << start_score << std::endl;
//	pose.dump_scored_pdb( "output/" + filename + "_pre_min.pdb", *score_fxn );

	pack::task::PackerTaskOP rotamerize_task( pack::task::TaskFactory::create_packer_task( pose ));
	rotamerize_task->initialize_from_command_line().restrict_to_repacking();
	rotamerize_task->set_bump_check( false );

	RotamerizeMover try_me_out( rotamerize_task );

	try_me_out.apply( pose );

	Real post_pack_interface_rmsd = core::scoring::rmsd_with_super( pose, start_pose, is_heavyatom );
	std::cout << "PostPack rmsd:  rmsd between pre and post min is " << post_pack_interface_rmsd << std::endl;
//	pose.dump_scored_pdb( "output/" + filename + "_rotamerized.pdb", *score_fxn );
	pose.dump_pdb( "output/" + filename + "_rotamerized.pdb" );

	return;
}

int
main( int argc, char * argv [] )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	devel::init( argc, argv );

	std::string pdb_path( option[ basic::options::OptionKeys::in::path::pdb ]()[1] );

	utility::vector1<file::FileName> file_names( option[ basic::options::OptionKeys::in::file::list ] );
	std::ifstream struct_list_file( file_names[1].name().c_str() );
	std::string pdb_code;

	struct_list_file >> pdb_code;
	while( !struct_list_file.eof() ){

		std::cout << "Working on file:  " << pdb_code << std::endl;

		pose::Pose target_pose;
		core::import_pose::pose_from_file( target_pose, pdb_path + pdb_code , core::import_pose::PDB_file);

		std::cout << "Total residue count is " << target_pose.size() << std::endl;

		rotamerize( pdb_code, target_pose );

		struct_list_file >> pdb_code;
	}


}


