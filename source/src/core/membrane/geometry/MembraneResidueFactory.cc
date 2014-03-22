// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/geometry/MembraneResidueFactory.cc
///
/// @brief 		Membrane Residue Factory - Creates residues of type MEM and EMB
/// @details 	Creates a membrane residue of type MEM with AA type virtual residue. Adds the residue by
///             either a specified jump or to the end of the fold tree and then makes the new residue the root.
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_geometry_MembraneResidueFactory_cc
#define INCLUDED_core_membrane_geometry_MembraneResidueFactory_cc

// Unit headers
#include <core/membrane/geometry/MembraneResidueFactory.hh>

// Project Headers
#include <core/conformation/membrane/definitions.hh>
#include <core/conformation/membrane/Exceptions.hh>

// Package Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>


// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>

// C++ Headers
#include <algorithm>
#include <string>
#include <cstdlib>
#include <cmath>

static basic::Tracer TR("core.membrane.geometry.MembraneResidueFactory");

namespace core {
namespace membrane {
namespace geometry {


/// @brief 	Default Constructor
MembraneResidueFactory::MembraneResidueFactory() {}

/// @brief  Default Destructor
MembraneResidueFactory::~MembraneResidueFactory() {}

/// @brief    Add a membrane definiiton residue to the pose
/// @details  Add a residue of type MEM to the pose at a Jump nres+1 and make the
///           mmebrane residue the root of the fold tree
void MembraneResidueFactory::add_membrane_residue(
													core::Vector& center,
													core::Vector& normal,
													core::Real& thickness,
													core::pose::Pose & pose,
													bool fullatom
																			) {
	using namespace core::membrane;
	using namespace core::conformation;
	using namespace core::chemical;

	TR << "Adding a membrane virtual residue to the pose" << std::endl;

	// Pose total residue
	int const nres = pose.total_residue();

	// Residue Type Info
	ResidueTypeSetCAP const & residue_set(
										  core::chemical::ChemicalManager::get_instance()->residue_type_set( fullatom ? core::chemical::FA_STANDARD : core::chemical::CENTROID )
										  );
	core::chemical::ResidueTypeCOPs const & rsd_type_list( residue_set->name3_map("MEM") );
	core::chemical::ResidueType const & membrane( *rsd_type_list[1] );

	// Create a new residue
	core::conformation::ResidueOP rsd( core::conformation::ResidueFactory::create_residue(membrane) );

	// Set Coordinates in the Residue Contianer
	core::Vector abs;
	abs.assign(0, thickness, 0);
	rsd->set_xyz(1, center);
	rsd->set_xyz(2, normal);
	rsd->set_xyz(3, abs);

	// Append residue by jump - no new chain
	pose.append_residue_by_jump( *rsd, 1, "", "", true ); // for lack of ft cycles...? changing to 1

	// Make the membrane residue the root of the fold tree
	int const vertex = nres+1;
	kinematics::FoldTree newF( pose.fold_tree() );
	newF.reorder(vertex);
	pose.fold_tree( newF );

	TR << "Membrane Residue added!" << std::endl;

}

/// @brief    Add a membrane definition residue to the pose
/// @details  Add a residue of type MEM to the pose at a Jump nres+1 and make the
///           membrane residue the root of the fold tree
void MembraneResidueFactory::add_embedding_residue(
		core::Vector& center,
		core::Vector& normal,
		core::Real& depth,
		core::pose::Pose & pose,
		Size jump,
		bool fullatom)
{
	using namespace core::membrane;
	using namespace core::conformation;
	using namespace core::chemical;

	TR << "Adding a membrane protein embedding virtual residue to the pose at jump anchor " << jump << std::endl;

	// Residue type information
	ResidueTypeSetCAP const & residue_set(core::chemical::ChemicalManager::get_instance()->residue_type_set( fullatom ?
			core::chemical::FA_STANDARD : core::chemical::CENTROID ));
	core::chemical::ResidueTypeCOPs const & rsd_type_list( residue_set->name3_map("EMB") );
	core::chemical::ResidueType const & embedding( *rsd_type_list[1] );

	// Create a new residue
	core::conformation::ResidueOP rsd( core::conformation::ResidueFactory::create_residue(embedding) );

	// Make a vector out of my depth
	// Note - depth vector must be created with x, z = 1. If x, z = 0, converting Cartesian to internal
	// coordinates in updating the atom tree will result in a divide by zero error for some coordinate sets.
	core::Vector abs;
	abs.assign(1, depth, 1);

	// Add my data to the virtual residue
	rsd->set_xyz(1, center);
	rsd->set_xyz(2, normal);
	rsd->set_xyz(3, abs);

	// Append residue to pose
	pose.append_residue_by_jump( *rsd, jump );

	// Reordering the foldtree to accommodate new residue
	kinematics::FoldTree newF( pose.fold_tree() );
	newF.reorder( 1 );
	pose.fold_tree( newF );

	TR << "Embedding residue added!" << std::endl;
}

} // geometry
} // membrane
} // core

#endif // INCLUDED_core_membrane_geometry_MembraneResidueFactory_cc

