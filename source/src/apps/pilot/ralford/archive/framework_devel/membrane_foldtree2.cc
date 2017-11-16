// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       apps/pilot/ralford/memrbane_foldtree.cc
///
/// @brief      Unit Test for Membrane Fold Tree Topologies
/// @details    Construct a membrane fodltree from a multi chain with an additional sequnce
///             of virtual residues strung to the end. This uses a function that mirrors a membrane
///             foldtree class included in my current branch rflaford12/membrane_frmwk_w_db
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (1/9/14)

// App Headers
#include <devel/init.hh>

// Project Headers
#include <core/kinematics/FoldTree.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>

#include <core/membrane/MembraneConformation.hh>

// Package Headers
#include <core/types.hh>

#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using basic::Error;
using basic::Warning;

static basic::Tracer TR( "apps.pilot.ralford.membrane_foldtree" );

/// @brief Load Membrane Pose
core::pose::PoseOP load_pose() {

	using namespace core::import_pose;
	using namespace core::pose;

	TR << "Loading 1afo from PDB" << std::endl;
	PoseOP pose = new Pose();
	pose_from_file( *pose, "test/core/membrane/io/1afo_test.pdb" , core::import_pose::PDB_file);

	return pose;
}

/// @brief Add Membrane and Embedding virtual residues
void add_vrts( core::pose::PoseOP pose ) {

	using namespace core::chemical;
	using namespace core::conformation;

	TR << "Adding membrane and embedding style virtual residues at the appropriate jumps" << std::endl;
	// Option Setting for residue type set
	ResidueTypeSetCAP const & residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ));
	core::chemical::ResidueTypeCOPs const & rsd_type_list1( residue_set->name3_map("VRT") );
	core::chemical::ResidueType const & virtuals( *rsd_type_list1[1] );

	// Create three virtual residues
	core::conformation::ResidueOP vrt1( core::conformation::ResidueFactory::create_residue(virtuals) );
	core::conformation::ResidueOP vrt2( core::conformation::ResidueFactory::create_residue(virtuals) );
	core::conformation::ResidueOP vrt3( core::conformation::ResidueFactory::create_residue(virtuals) );

	// Setting up traditional membrane-style coordinates (vrt1 simulates the membrane)
	core::Vector center(0, 0, 0);
	core::Vector normal(0, 0, 1);
	core::Vector depth( 1, 30.0, 1);

	vrt1->set_xyz( 1, center );
	vrt1->set_xyz( 2, normal );
	vrt1->set_xyz( 3, depth );

	// doing this for 2 and 3 as a quick fix
	vrt2->set_xyz( 3, depth );
	vrt3->set_xyz( 3, depth );

	// Add virtual residues to the existing pose
	pose->append_residue_by_jump( *vrt1, 1, "", "", true);
	pose->append_residue_by_jump( *vrt2, 1 );
	pose->append_residue_by_jump( *vrt3, 41 );

	// Add an nterm patch

	return;
}

/// @brief Modify the current fold tree topology
void modify_foldtree( core::pose::PoseOP pose ) {

	using namespace core::membrane;

	TR << "Constructing a new conformation" << std::endl;

	// Setup the embedding residue map
	std::map< int, int > embres_map;
	embres_map.insert( std::pair< int, int >( 1, 82 ) );
	embres_map.insert( std::pair< int, int >( 2, 83 ) );

	// Setup the membrane root
	int membrane = 81;

	MembraneConformationOP mp_conf( new MembraneConformation( pose->conformation(), embres_map, membrane ));
}

/// @brief Main Function
int main( int argc, char* argv[] )
{
	try {

		// Initialize Options System, RG, and All Factory_Registrators
		devel::init(argc, argv);

		TR << "Pilot App: Membrane Fold Tree" << std::endl;
		TR << "Author: Rebecca Alford" << std::endl;
		TR << "Testing membrane foldtree topologies" << std::endl;

		// Set up a pose from pdb
		core::pose::PoseOP pose = load_pose();

		// Add virtual residues
		add_vrts(pose);

		// Integrate and print consecutive fold trees
		modify_foldtree(pose);

		TR << "Done!" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
}

