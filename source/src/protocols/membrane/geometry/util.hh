// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		protocols/membrane/geometry/util.hh
///
/// @brief 		Utility methods for membrane framework
/// @details 	Utility methods include determining center of mass (moved down in the tree)
///				and adjusting normal parameters for visualization.
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_geometry_util_hh
#define INCLUDED_protocols_membrane_geometry_util_hh

// Package Headers
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/membrane/MembraneInfo.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>

#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

// C++ Headers
#include <cstdlib>
#include <cmath>

namespace protocols {
namespace membrane {
namespace geometry {

using namespace core; 

//////////////// Utility Functions from Docking Protocol - Geometry Util for Center of Mass ////////////////

/// @brief      Center of Mass
/// @details    Calculates the center of mass of a pose - Stop and start positions (or residues)
///             used ot find the starting and finishing locations
///				the start and stop positions (or residues) within the pose are used to
///				find the starting and finishing locations
///
/// @author     Monica Berrondo, Modified by Javier Castellanos and Rebecca Alford
numeric::xyzVector< core::Real>
center_of_mass(
			   core::pose::Pose const & pose,
			   core::SSize const start,
			   core::SSize const stop
			   );

/// @brief      Residue Center of Mass
/// @details    Calcualte the center of mass of a pose.
///
/// @author     Monica Berrondo, Modified by Javier Castellanos and Rebecca Alford
core::SSize
residue_center_of_mass(
					   core::pose::Pose const & pose,
					   core::SSize const start,
					   core::SSize const stop
					   );


/// @brief      Return nearest residue
/// @details    Find the residue nearest some position passed in (normally a center of mass)
///
/// @author     Monica Berrondo, Modified by Javier Castellanos and Rebecca Alford
core::SSize
return_nearest_residue(
					   core::pose::Pose const & pose,
					   core::SSize const begin,
					   core::SSize const end,
					   core::Vector center
					   );
				
/// @brief		Get z-coord and chainID
/// @details	Helper function that creates input for SpanningTopology
///				which is not built at the time the Pose is built
///				returns a pair of vectors:
///				vector1 is z-coord of CA atoms of the pose
///				vector2 is chainID of CA atoms of the pose
std::pair< utility::vector1< Real >, utility::vector1< Real > > get_chain_and_z( pose::PoseOP pose );

/// @brief Normalize normal vector to length 15 for visualization
void membrane_normal_to_length_15( pose::Pose & pose );


} // geometry
} // membrane
} // protocols

#endif // INCLUDED__membrane_geometry_util_hh

