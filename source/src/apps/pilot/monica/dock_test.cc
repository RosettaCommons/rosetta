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
/// @author Modified by Sergey Lyskov


// libRosetta headers
#include <protocols/docking/DockingProtocol.hh>

#include <numeric/random/random.hh>


#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>


#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>


#include <core/types.hh>


#include <protocols/rigid/RB_geometry.hh>

#include <core/io/pdb/build_pose_as_is.hh>

#include <core/io/pdb/pdb_writer.hh>

//#include <core/pose/PScene.hh>

#include <devel/init.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>


//silly using/typedef

using namespace core;

using utility::vector1;


typedef std::map< std::string, std::map< std::string, numeric::xyzVector< Real > > > Coords;

typedef vector1< std::string > Strings;
///////////////////////////////////////////////////////////////////////////////
// some silly helper routines:


///////////////////////////////////////////////////////////////////////////////
std::string readFile(std::string fname)
{
	Size fsize;
	std::string res;

	std::ifstream file(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
	if( file.is_open() )  {
		fsize = file.tellg();
		res.resize( fsize );
		file.seekg(0, std::ios::beg);
		file.read(&res[0], fsize);
		file.close();
	}
	else std::cout << "file not found!";
	return res;
}

///////////////////////////////////////////////////////////////////////////////
void
rb_test ()
{
	using namespace pose;
	using namespace kinematics;
	using namespace protocols::moves;
	using namespace protocols::geometry;

	using core::import_pose::pose_from_file;

	// PHIL temporarily hacking for Jeff until file_data is rehabilitated

	//chemical::ResidueTypeSetCAP residue_set
	//	( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );

	std::cout << "---------------------"<< std::endl;

	std::cout << "\n\n\n\n\n-------------------------------------\n";
	// TODO get rid of FileData
	//std::string data1 = readFile("test_in.pdb");
	std::string data1 = readFile("dock_protein1.pdb");
	std::string data2 = readFile("dock_protein2.pdb");

	core::import_pose::FileData fd_one = core::io::pdb::PDB_DReader::createFileData(data1);
	core::import_pose::FileData fd_two = core::io::pdb::PDB_DReader::createFileData(data2);

	core::pose::Pose pose;
	// somehow create a pose from the two pdbs that were read in
	//fd_one.build_pose(pose, residue_set);
	core::import_pose::pose_from_file( pose, "dock_protein1.pdb", core::import_pose::PDB_file);
	std::cout << "total_residue = " << pose.total_residue() << "\n";
	std::cout << "###########################################" << std::endl;

	core::pose::Pose tmp_pose;
	//fd_two.build_pose(tmp_pose, residue_set);
	core::import_pose::pose_from_file( tmp_pose, "dock_protein2.pdb", core::import_pose::PDB_file);

	std::cout << "total_residue = " << pose.total_residue() << "\n";
	int cutpoint ( pose.total_residue() );
	pose.dump_pdb( "pose1.pdb" );
//	pose.copy_segment( tmp_pose.total_residue(), tmp_pose, 1, pose.total_residue()+1 );
//	TO DO fix me!!
	for ( Size i=1; i<=tmp_pose.total_residue(); ++i ) {
		conformation::ResidueCOP new_rsd = tmp_pose.residue(i).clone();
		if ( i == 1 ) {
			// since the first residue is a terminus variant, it doesn't know how to connect to a preceding residue
			// ergo, attach by a jump
			pose.append_residue_by_jump( *new_rsd, cutpoint );
		} else {
			pose.append_residue_by_bond( *new_rsd );
		}
	}

	std::cout << "total_residue2 = " << pose.total_residue() << "\n";
	pose.dump_pdb( "pose2.pdb" );

	int const nres( pose.total_residue() );
	int dock_jump ( 0 );
	assert( nres > 40 );
	// this actually needs to be a one_jump_tree
	// one_jump_tree exists in old rosetta, but not in minirosetta
	// we need to see whether to create a one_jump_tree or just create
	// the tree using tree_from_jumps_and_cuts
	{ // set up the foldtree
		FoldTree f( nres );
		int const jump_pos1 ( core::pose::residue_center_of_mass( pose, 1, cutpoint ) );
		int const jump_pos2 ( core::pose::residue_center_of_mass( pose, cutpoint+1, pose.total_residue() ) );
		dock_jump = f.new_jump( jump_pos1, jump_pos2, cutpoint );
		// set this jump as the docking jump
		f.reorder( jump_pos1 );
		pose.fold_tree( f );
	}
	std::cout << " jump tree created" << std::endl;

	// run docking
	protocols::docking::DockingProtocol dock( dock_jump );
	dock.apply( pose );

	pose.dump_pdb( "finalstructure.pdb" );

	// store and output
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief main docking protocol
///
/// @details
/// Main function for docking. Includes the following steps:
///      0) prepack mode: prepare a starting structure for later runs
///   OR:
///      1) perturbation of decoy (see docking_perturb_decoy): changes
///         orientation of docking partners
///      2) low-resolution search:
///         refine structure with rigid-body, centroid-mode MC cycles
///      3) high-resolution search:
///         further refinement of structure in fullatom mode

int
main( int argc, char * argv [] )
{
	try{
	using namespace core;

	// initialize core
	devel::init(argc, argv);

// 	rb_test( residue_set );
	rb_test();

	std::cout << "Done! -------------------------------\n";
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

