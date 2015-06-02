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
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/membrane/MembraneInfo.fwd.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <protocols/membrane/geometry/EmbeddingDef.fwd.hh>

#include <core/pose/Pose.fwd.hh>
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
using namespace core::conformation::membrane;

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
std::pair< utility::vector1< Real >, utility::vector1< Real > > get_chain_and_z( pose::Pose const & pose );

/// @brief		Get secondary structure of the pose
/// @details	Helper function that gets the secondary structure vector from the pose
utility::vector1< char > get_secstruct( pose::Pose & pose );


/// @brief Compute Membrane Center/Normal from Membrane Spanning
/// topology
void compute_structure_based_embedding(
			  pose::Pose const & pose,
			  SpanningTopology const & topology,
			  Vector & center,
			  Vector & normal
			  );

/// @brief Compute Membrane Center/Normal from Membrane Spanning
/// topology, uses topology from MembraneInfo
void compute_structure_based_embedding(
			   pose::Pose const & pose,
			   Vector & center,
			   Vector & normal
			   );

/// @brief Compute Membrane Center/Normal from Membrane Spanning
/// topology
protocols::membrane::geometry::EmbeddingDefOP
compute_structure_based_embedding( pose::Pose const & pose, SpanningTopology const & topo );

/// @brief Compute Membrane Center/Normal from Membrane Spanning
/// topology, uses topology from MembraneInfo
protocols::membrane::geometry::EmbeddingDefOP
compute_structure_based_embedding( pose::Pose const & pose );

/// @brief Check reasonable range of vector
void check_vector( core::Vector const vector );

/// @brief Average EmbeddingDefs as they are (without vector inversion accounting for topology)
/// @details Get average center and normal from a vector of EmbeddingDefs
EmbeddingDefOP average_embeddings( utility::vector1< EmbeddingDefOP > const parts );

/// @brief Average EmbeddingDefs after first inverting some vectors accounting for topology
/// @details Get average center and normal from a vector of EmbeddingDefs
EmbeddingDefOP average_antiparallel_embeddings( utility::vector1< EmbeddingDefOP > const parts );

/// @brief Normalize normal vector to length 15 for visualization
void membrane_normal_to_length_15( pose::Pose & pose );

/// @brief Set membrane residue to root of foldtree
/// @details Requires MembraneInfo to be constructed beforehand;
///			 use AddMembraneMover to do that
void reorder_membrane_foldtree( pose::Pose & pose );

/// @brief Calculates translation axis lying in the membrane (= projection axis
///			between embedding centers)
core::Vector const membrane_axis( pose::Pose & pose, int jumpnum );

/// @brief Splits the SpanningTopology object into two objects, depending on
///			given jump number
/// @details This is useful for calculating an embedding for parts of the
///			structure: this can now easily be accomplished by creating two
///			empty topology objects, call this function, and then use both topology
///			objects and subposes to call compute_structure_based_membrane_embedding
///			BEWARE: this does not work for splitting topology by spans! It only works
///			chainwise
void split_topology_by_jump(
		pose::Pose const & pose,		// full pose
		Size const jumpnum,				// jump number to split on
		SpanningTopology const & topo,	// topology to split
		pose::Pose & pose_up,			// upstream partner after pose splitting
		pose::Pose & pose_down,			// downstream partner after pose splitting
		SpanningTopology & topo_up,		// topology of upstream pose
		SpanningTopology & topo_down	// topology of downstream pose
		);

/// @brief Splits the SpanningTopology object into two objects, depending on
///			given jump number
/// @details This doesn't shift the topology to start at 1 for each partner, it
///				remains exactly the same as it would be for the complete pose, just split
///			BEWARE: this does not work for splitting topology by spans! It only works
///			chainwise
void split_topology_by_jump_noshift(
		pose::Pose const & pose,		// full pose
		Size const jumpnum,				// jump number to split on
		SpanningTopologyOP topo,	// topology to split
		SpanningTopologyOP topo_up,		// topology of upstream pose
		SpanningTopologyOP topo_down	// topology of downstream pose
		);

/// @brief Update embedding of the partners after a move
/// @details Requires the jump number between the partners, the topology will
///				be taken from MembraneInfo and will be split accordingly; up and
///				down means upstream and downstream
void update_partner_embeddings( pose::Pose const & pose, Size const jumpnum, EmbeddingDef & emb_up, EmbeddingDef & emb_down );

} // geometry
} // membrane
} // protocols

#endif // INCLUDED__membrane_geometry_util_hh

