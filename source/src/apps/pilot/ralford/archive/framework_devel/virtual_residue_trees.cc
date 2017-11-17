// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       apps/pilot/ralford/membrane_integration.cc
///
/// @brief      Top-Level Unit Test for the Membrane Protein Factory
/// @details    The purpose of this application is to test the membrane protein factory
///             initialization code including all external dependencies which cannot be tested
///             in JD2, Resource Manager, and the Pose cache. This can also serve as the integration
///             test for memrane protein initialization.
///
/// @note       This test is highly coupled to the xml file membrane.xml
/// @note       This test contains _no_ scoring integration!!!
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (1/2/14)

// App Headers
#include <devel/init.hh>

// Project Headers
#include <core/kinematics/FoldTree.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

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

static basic::Tracer TR( "apps.pilot.ralford.virtual_residue_tree" );

/// @brief Build Toy Pose
core::pose::PoseOP build_toy_pose() {

	using namespace core::conformation;
	using namespace core::chemical;
	using namespace core::pose;
	using namespace core::kinematics;

	// Make a new pose
	PoseOP pose = new Pose();

	// Option Setting for residue type set
	ResidueTypeSetCAP const & residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ));

	// Set up options for adding alanines
	core::chemical::ResidueTypeCOPs const & rsd_type( residue_set->name3_map("ALA") );
	core::chemical::ResidueType alanine = *rsd_type[1];

	// Make a bunch of alanines
	core::conformation::ResidueOP rsd1( core::conformation::ResidueFactory::create_residue(alanine) );
	core::conformation::ResidueOP rsd2( core::conformation::ResidueFactory::create_residue(alanine) );
	core::conformation::ResidueOP rsd3( core::conformation::ResidueFactory::create_residue(alanine) );
	core::conformation::ResidueOP rsd4( core::conformation::ResidueFactory::create_residue(alanine) );

	// Append these residues first to make a small pose
	pose->append_residue_by_bond( *rsd1 );
	pose->append_residue_by_bond( *rsd2 );
	pose->append_residue_by_bond( *rsd3 );
	pose->append_residue_by_bond( *rsd4 );

	// Set up options for adding a string of virtual residues afterwards
	core::chemical::ResidueTypeCOP rsd_type( residue_set->get_representative_type_name3("VRT") );
	core::chemical::ResidueType virtuals = *rsd_type;

	// Make a few virtual residues
	core::conformation::ResidueOP rsd5( core::conformation::ResidueFactory::create_residue(virtuals) );
	core::conformation::ResidueOP rsd6( core::conformation::ResidueFactory::create_residue(virtuals) );
	core::conformation::ResidueOP rsd7( core::conformation::ResidueFactory::create_residue(virtuals) );
	core::conformation::ResidueOP rsd8( core::conformation::ResidueFactory::create_residue(virtuals) );

	// Add a string of virtual residues
	pose->append_residue_by_jump( *rsd5, 1, "", "", true );
	pose->append_residue_by_jump( *rsd6, 2, "", "", false );
	// will add the rest later

	TR << "Printing information about my newly constructed pose!" << std::endl;
	TR << "Number of residues: " << pose->size() << std::endl;
	TR << "Number of chains: " << pose->conformation().num_chains() << std::endl;
	TR << "Printing the resulting foldtree" << std::endl;
	pose->fold_tree().show(std::cout);

	TR << "Modifying the foldtree topology from scratch to have a membrane-like topology" << std::endl;
	// Make a new fold tree with the desired topology
	FoldTreeOP ft = new FoldTree();
	ft->add_edge( 1, 5, 1);
	ft->add_edge( 1, 4, -1);
	ft->add_edge( 1, 6, 2);

	// Try reorder
	ft->reorder( 5 );

	// Add variant types!!!
	core::pose::add_variant_type_to_pose_residue( *pose, "LOWER_TERMINUS", 1);
	core::pose::add_variant_type_to_pose_residue( *pose, "UPPER_TERMINUS", 4);

	// Take the fold tree and try to reset it
	pose->fold_tree( *ft );

	// Reorder


	TR << "Showing the new foldtree" << std::endl;
	pose->fold_tree().show(std::cout);

	return pose;

}

/// @brief Main
int main( int argc, char* argv[] )
{

	try {

		// Initialize Options System, RG, and All Factory_Registrators
		devel::init(argc, argv);

		TR << "Building a toy pose which contains a string of alanine residues and then a proceeding string of virtual residues" << std::endl;
		build_toy_pose();

		TR << "Done!" << std::endl;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

