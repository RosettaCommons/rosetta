// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/util.hh
///
/// @brief      Utility methods for working with proteins in the membrane
/// @detiails   Several groups of utilities for working in the memrbane environment
///                 * Calculate RMSD between the transmembrane domains of two poses
///                   (with or without superimposition)
///                 * Calculate the tilt of the protien relative to the membrane
///                   normal
///                 * Safety checks and convenience methods for working with
///                   membrane foldtrees
///                 * Utility for accessing DSSP secstruc and z coordinates
///                 * Calculate protein embedding based on the structure
///                 * Split topology by jump, and other multi-chain (or partner)
///                   functions
///                 * Methods for reducing the amount of code duplication in RosettaScripts
///                   and init_from_cmd functions in the membrane framework
///
/// NOTE: All of these methods require a RosettaMP framework pose or eventually
/// may require this. Use pose.conformation().is_membrane() for safety checks!
///
/// Last Modified: 7/9/15
/// @author Rebecca faye Alford (rfalford12@gmail.com)
/// @author JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_util_hh
#define INCLUDED_protocols_membrane_util_hh

// Package Headers
#include <core/conformation/membrane/MembraneInfo.fwd.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <protocols/membrane/geometry/EmbeddingDef.fwd.hh>
#include <protocols/membrane/geometry/Embedding.fwd.hh>

// Project Headers
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.fwd.hh>

#include <core/kinematics/FoldTree.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/tag/Tag.hh>

#include <utility/vector1.hh>

#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

// C++ headers
#include <cstdlib>
#include <cmath>

namespace protocols {
namespace membrane {

/////////////////////////////////////////////////////////////////////////
// Methods for calculating rmsds between protein transmembrane regions //
/////////////////////////////////////////////////////////////////////////

/// @brief Compute backbone RMSD between TM regions - don't superimpose
/// @details Calculate the rmsd between backbone atoms (N, CB, CA, O)
/// in the transmembrane regions, as defined by the spanning topology
/// object Do not superimpose the poses. Takes a native pose and
/// current pose
core::Real
mem_bb_rmsd_no_super(
	core::pose::Pose & native_pose,
	core::pose::Pose & pose
);

/// @brief Compute all-atom RMSD between TM regions - don't superimpose
/// @details Calculate the rmsd between all atoms in the pose in the
/// transmembrane regions, as defined by the spanning topology object.
/// Do not superimpose the poses. Takes a native pose & current pose
core::Real
mem_all_atom_rmsd_no_super(
	core::pose::Pose & native_pose,
	core::pose::Pose & pose
);

/// @brief Compute backbone RMSD between TM regions - do superimpose
/// @details Calculate the rmsd between backbone atoms (N, CB, CA, O)
/// in the transmembrane regions, as defined by the spanning
/// topology object Superimpose the poses. Takes a native pose and
/// current pose
core::Real
mem_bb_rmsd_with_super(
	core::pose::Pose & native_pose,
	core::pose::Pose & pose
);

/// @brief Compute all-atom RMSD between TM regions - do superimpose
/// @details Calculate the rmsd between all atoms in the pose in the
/// transmembrane regions, as defined by the spanning topology object.
/// Superimpose the poses. Takes a native pose & current pose
core::Real
mem_all_atom_rmsd_with_super(
	core::pose::Pose & native_pose,
	core::pose::Pose & pose
);

//////////////////////////////////////////////////////////////////
// Methods for calculating tilt of helices relative to membrane //
//////////////////////////////////////////////////////////////////

/// @brief Calculate tilt of a TM span relative to the membrane normal
/// @details Given a transmembrane span #, calculate the angle between the
/// axis through the helix and the membrane normal. Works for relatively
/// straight helices but less accurate for kinks. Takes a pose & span number.
core::Real
calc_helix_tilt_angle( core::pose::Pose & pose, core::Size span_no );

/// @brief Determine the axis used to define a single TM Helix
/// @details Using the COM of the helix start & end position, calculate a helix
/// describing its geometry relative to the memrbane normal. Takes a pose &
/// span number. Not a good approx for helices with kinks.
/// TODO: CODE INSIDE SHOULD BE REPLACED WITH SPAN EMBEDDING CALCULATIONS!
///       THAT'S WHAT THEY ARE THERE FOR! (JKLeman)
core::Vector
calc_helix_axis( core::pose::Pose & pose, core::Size span_no );

/// @brief Calculate center of mass between 3 xyz coords
/// @details Given three xyz vectors, calculate the center of mass
/// and return a vector. Helper method to calc_helix axis.
core::Vector
com( core::Vector a, core::Vector b, core::Vector c );

/// @brief Calculate the RMSD between a helix tilt angle & reference
/// @details Given a reference angle and measured angle, calculate the
/// root mean square deviation between the two single values. Takes
/// the measured tilt angle and reference angle (typically from experiment)
core::Real
calc_angle_rmsd( core::Real measured_angle, core::Real ref_angle );

/// @brief Calculate tilt angle and distance from membrane center
/// @details Computes the tilt angle and distance from membrane center for the
///   protein embedding center and normal
utility::vector1< core::Real > pose_tilt_angle_and_center_distance( core::pose::Pose & pose );

////////////////////////////////////////////////////////////////
// Safety checks & convenience methods for membrane foldtrees //
////////////////////////////////////////////////////////////////

/// @brief Determine whether the membrane is modeled as fixed
/// @details Based on the setup of the foldtree, determined whether
/// the membrane is currently fixed, meaning it is setup as the
/// root in the FoldTree and has no upstream children. Takes a pose.
bool
is_membrane_fixed( core::pose::Pose & pose );

/// @brief Determine whether membrane can move on its own
/// @details Based on the setup of the FoldTree, determine whether
/// the membrane is moveable, but when moved, won't cause anything in the
/// protein to move (i.e. independently moveable). Takes a pose.
bool
is_membrane_moveable_by_itself( core::pose::Pose & pose );

/// @brief Set membrane residue to root of foldtree
/// @details Naively sets the root of the foldtree to be the membrane
/// residue. Should perform checks before doing this!
void reorder_membrane_foldtree( core::pose::Pose & pose );

/// @brief Create a membrane foldtree with an interface
/// @details Currently only works for two-body-docking. Both partners can have
///   multiple chains in any order, the anchoring happens at the TM COM
///   of each chain
///
///       __________________________________________
///      |________________  _________________       |
///      |________   iJ   ||________         |      |
///      |        |       ||        |        |      |
/// -------  -------  -------  -------  -------  M=root
///  chain1   chain2   chain3   chain4 ...
///
/// ^--------------^  ^-----------------------^
///     partner 1             partner 2
///
///  iJ = interface jump, will be returned from the function
///
core::Size create_membrane_docking_foldtree_from_partners( core::pose::Pose & pose, std::string const partners );

/// @brief Create membrane foldtree from scratch
/// @details The foldtree is setup such that the membrane is at the root and
///   anchored at the first chain COM residue with jumps from the
///   first chain COM to each chain COM; requires the membrane to be present
///   Returns the root anchoring point, i.e. rsd nearest chain COM of chain 1
///
///     ________________________________
///    |__________________________      |
///    |_________________         |     |
///    |________         |        |     |
///    |        |        |        |     |
/// -------  -------  -------  -------  M=root
///  chain1   chain2   chain3   chain4 ...
core::Size create_membrane_foldtree_anchor_com( core::pose::Pose & pose );

/// @brief Create membrane foldtree from scratch
/// @details The foldtree is setup such that the membrane is at the root and
///   anchored at the first chain TRANSMEMBRANE COM residue with jumps from the
///   first chain COM to each chain TRANSMEMBRANE COM;
///   requires the membrane to be present
///   Returns the root anchoring point, i.e. rsd nearest chain TM COM of
///   chain 1
///
///       ________________________________
///      |__________________________      |
///      |_________________         |     |
///      |________         |        |     |
///      |        |        |        |     |
/// -------  -------  -------  -------  M=root
///  chain1   chain2   chain3   chain4 ...
core::Size create_membrane_foldtree_anchor_tmcom( core::pose::Pose & pose );

/// @brief Create membrane foldtree from scratch
/// @details The foldtree is setup such that the membrane is at the root and
///   anchored at the residue closest to the pose TM COM with jumps from there
///   to each chain TRANSMEMBRANE COM; requires the membrane to be present;
///   Returns the root anchoring point, i.e. rsd nearest pose TM COM that
///   can be in any chain
///
///                _______________________
///               |_________________      |
///               |________         |     |
///       ________|        |        |     |
///      |        |        |        |     |
/// -------  -------  -------  -------  M=root
///  chain1   chain2   chain3   chain4 ...
core::Size create_membrane_foldtree_anchor_pose_tmcom( core::pose::Pose & pose );

/// @brief Create foldtree for multiple docking partners
/// @details The foldtree is setup such that the membrane is at the root and
///   anchored at the residue closest to the TM COM of the first chain.
///   The partner setup defines how the foldtree will look like: partners
///   are anchored at the first chain, and multiple chains in a partner
///   are anchored to the first chain of that partner, example below.
///   Returns all interface jump numbers.
///   Example: A_BCD_E
///
///      __________________________________________
///     |________________________________          |
///     |          _______________       |         |
///     |_______  |_______        |      |         |
///     |       | |       |       |      |         |
/// -------  -------  -------  -------  -------  M=root
///  chain1   chain2   chain3   chain4   chain5...
///
/// ^-----^  ^-----------------------^  ^-----^
/// partner1          partner2          partner3
///
utility::vector1< core::Size > create_membrane_multi_partner_foldtree_anchor_tmcom( core::pose::Pose & pose, std::string partner );

/// @brief Helper function to create membrane foldtrees
/// @details The anchors vector is a vector of anchor residues in all chains,
///   one per chain. This function assumes that the first entry in the
///   vector is the root anchor point to which all other chains are
///   connected;
///       ________________________________
///      |__________________________      |
///      |_________________         |     |
///      |________         |        |     |
///      |        |        |        |     |
/// -------  -------  -------  -------  M=root
///  chain1   chain2   chain3   chain4 ...
///
void create_membrane_foldtree_from_anchors( core::pose::Pose & pose, utility::vector1< core::Size > anchors );

/// @brief Helper function to create a specific membrane foldtree
/// @details I am hijacking xyzVectors to hold the jumps that need to be
///   created in the foldtree: xyz = ( rsd1, rsd2, cutpoint )
///   THE JUMP NUMBERS WILL BE CONSECUTIVE, ACCORDING TO THE VECTOR1
core::Size create_specific_membrane_foldtree( core::pose::Pose & pose, utility::vector1< core::Vector > anchors );

/// @brief Helper function to create membrane foldtrees
/// @details Returns the residues closest to the COMs for each chain
utility::vector1< core::Size > get_anchor_points_for_tmcom( core::pose::Pose & pose );

/// @brief Setup foldtree from scratch
/// @details The foldtree is setup such that the residue closest to the
///   COM is at the root, with jumps from there to each chain COM;
///   requires the membrane to be present;
///   Returns the root anchoring point, i.e. rsd nearest pose COM that
///   can be in any chain
///                _________________
///               |________         |
///       ________|        |        |
///      |        |        |        |
/// -------  -------  -------  -------
///  chain1   chain2   chain3   chain4 ...
///               ^
///       root anchor point
///
core::Size setup_foldtree_pose_com( core::pose::Pose & pose );

/// @brief Helper function to setup foldtrees
/// @details The anchors vector is a vector of anchor residues in all chains,
///   one per chain. This function assumes that the first entry in the
///   vector is the root anchor point to which all other chains are
///   connected;
///       __________________________
///      |_________________         |
///      |________         |        |
///      |        |        |        |
/// -------  -------  -------  -------
///  chain1   chain2   chain3   chain4 ...
///      ^
/// root anchor point
///
void setup_foldtree_from_anchors( core::pose::Pose & pose, utility::vector1< core::Size > anchors );

///////////////////////////////////////////////////////////
// Utilities for accessing dssp, z coords and chain info //
///////////////////////////////////////////////////////////

/// @brief Grab the z-coordinates and chainIDs from the entire pose
/// @details From the pose, grab all of the z_coords of CA atoms and
/// chain IDs, currently used for spanning topology construction.
/// Returns a std::pair of two vectors: the first a vector1 of z
/// coordinates and the second a vector1 of chainIDs for CA atoms
std::pair< utility::vector1< core::Real >, utility::vector1< core::Real > > get_chain_and_z( core::pose::Pose const & pose );

/// @brief  Get dssp defined secondary structure from the pose
/// @details Given a pose, grab a vector of characters describing the secondary
/// structure at each residue position in the pose, defined by DSSP
utility::vector1< char > get_secstruct( core::pose::Pose & pose );


///////////////////////////////////////////////////////////////////
// Methods for calculating the protein embedding in the membrane //
///////////////////////////////////////////////////////////////////

/// @brief Compute Membrane Center/Normal from Membrane Spanning
/// topology
void
compute_structure_based_embedding(
	core::pose::Pose const & pose,
	core::conformation::membrane::SpanningTopology const & topology,
	core::Vector & center,
	core::Vector & normal
);

/// @brief Compute Membrane Center/Normal from Membrane Spanning
/// topology, uses topology from MembraneInfo
void
compute_structure_based_embedding(
	core::pose::Pose const & pose,
	core::Vector & center,
	core::Vector & normal
);

/// @brief Compute Membrane Center/Normal from Membrane Spanning
/// topology
protocols::membrane::geometry::EmbeddingDefOP
compute_structure_based_embedding(
	core::pose::Pose const & pose,
	core::conformation::membrane::SpanningTopology const & topo
);

/// @brief Compute Membrane Center/Normal from Membrane Spanning
/// topology, uses topology from MembraneInfo
protocols::membrane::geometry::EmbeddingDefOP
compute_structure_based_embedding( core::pose::Pose const & pose );

/// @brief Compute embeddings by chain
/// @details The embeddings can be computed either from pose and topology or they
///   can be optimized differently; the function correlates each EmbeddingDef
///   object in embeddings with a span object in the pose's topology;
///   The bool means whether embeddings in Embedding object are
///   antiparallel or not
protocols::membrane::geometry::EmbeddingOP
compute_embeddings_by_chain( core::pose::Pose const & pose );

/// @brief Average EmbeddingDefs as they are (without vector inversion accounting for topology)
/// @details Get average center and normal from a vector of EmbeddingDefs
protocols::membrane::geometry::EmbeddingDefOP
average_embeddings( utility::vector1< protocols::membrane::geometry::EmbeddingDefOP > const parts );

/// @brief Average EmbeddingDefs after first inverting some vectors accounting for topology
/// @details Get average center and normal from a vector of EmbeddingDefs
protocols::membrane::geometry::EmbeddingDefOP
average_antiparallel_embeddings(
	utility::vector1< protocols::membrane::geometry::EmbeddingDefOP > const parts );

/// @brief Update embedding of the partners after a move
/// @details Requires the jump number between the partners, the topology will
///    be taken from MembraneInfo and will be split accordingly; up and
///    down means upstream and downstream
void
update_partner_embeddings(
	core::pose::Pose const & pose,
	core::Size const jumpnum,
	protocols::membrane::geometry::EmbeddingDef & emb_up,
	protocols::membrane::geometry::EmbeddingDef & emb_down );

/// @brief Pose transmembrane center-of-mass
/// @details Gets the coordinates of the TM span center-of-mass
///   This only looks at the span start and end residues for
///   calculation of the TM span COM, this should be faster than the
///   real thing though
core::Vector pose_tm_com( core::pose::Pose const & pose );

/// @brief Chain center-of-mass
/// @details Gets the coordinates of the chain center-of-mass
core::Vector chain_com( core::pose::Pose const & pose, core::Size chain );

/// @brief Chain center-of-mass of TM regions
/// @details Gets the coordinates of the chain center-of-mass but only the TM regions
core::Vector chain_tm_com( core::pose::Pose const & pose, core::Size chain );

/// @brief Residue closest to pose transmembrane center-of-mass
/// @details Gets the coordinates of the residue closest to TM span center-of-mass
///   This only looks at the span start and end residues for
///   calculation of the TM span COM, this should be faster than the
///   real thing though
core::Size rsd_closest_to_pose_tm_com( core::pose::Pose const & pose );

/// @brief Residue closest to chain center-of-mass
/// @details Gets the residue number closest to the chain center-of-mass
core::Size rsd_closest_to_chain_com( core::pose::Pose const & pose, core::Size chain );

/// @brief Residue closest to chain TM center-of-mass
/// @details Gets the residue number closest to the chain TM center-of-mass
core::Size rsd_closest_to_chain_tm_com( core::pose::Pose const & pose, core::Size chainid );


//////////////////////
// Membrane vectors //
//////////////////////

/// @brief Check reasonable range of vector
void check_vector( core::Vector const vector );

/// @brief Normalize normal vector to length 15 for visualization
void membrane_normal_to_length_15( core::pose::Pose & pose );

/// @brief Calculates translation axis lying in the membrane (= projection axis
///   between embedding centers)
core::Vector const membrane_axis( core::pose::Pose & pose, int jumpnum );

//////////////////////////////////////////////////////////////
// Methods for working with multiple partners and/or chains //
//////////////////////////////////////////////////////////////

/// @brief Splits the SpanningTopology object into two objects, depending on
/// given jump number
/// @details This is useful for calculating an embedding for parts of the
/// structure: this can now easily be accomplished by creating two empty topology
/// objects, call this function, and then use both topology objects and subposes
/// to call compute_structure_based_membrane_embedding
/// BEWARE: this does not work for splitting topology by spans! It only works
/// chainwise
void split_topology_by_jump(
	core::pose::Pose const & pose,  // full pose
	core::Size const jumpnum,    // jump number to split on
	core::conformation::membrane::SpanningTopology const & topo, // topology to split
	core::pose::Pose & pose_up,   // upstream partner after pose splitting
	core::pose::Pose & pose_down,   // downstream partner after pose splitting
	core::conformation::membrane::SpanningTopology & topo_up,  // topology of upstream pose
	core::conformation::membrane::SpanningTopology & topo_down // topology of downstream pose
);

/// @brief Splits the SpanningTopology object into two objects, depending on
/// given jump number
/// @details This doesn't shift the topology to start at 1 for each partner, it
/// remains exactly the same as it would be for the complete pose, just split
/// BEWARE: this does not work for splitting topology by spans! It only works
/// chainwise
void split_topology_by_jump_noshift(
	core::pose::Pose const & pose,  // full pose
	core::Size const jumpnum,    // jump number to split on
	core::conformation::membrane::SpanningTopologyOP topo,        // topology to split
	core::conformation::membrane::SpanningTopologyOP topo_up,  // topology of upstream pose
	core::conformation::membrane::SpanningTopologyOP topo_down // topology of downstream pose
);

/// @brief Split topology by chain
/// @details Split topology by chain and give vector of topology objects
utility::vector1< core::conformation::membrane::SpanningTopologyOP >
split_topology_by_chain_noshift(
	core::pose::Pose const & pose,
	core::conformation::membrane::SpanningTopologyOP const topo
);


/////////////////////////////////////////////
// Methods for reading center/normal in IO //
/////////////////////////////////////////////

/// @brief Read in a user provided center/normal pair from RosettaScripts
/// @details Given an XML tag from a RosettaScript read in a center & normal
/// option into two xyzVector objects. This method is intended to reduce duplication
/// accross membrane framework movers that use the same tricks for vectors.
/// Takes two Vector references and a Tag&
void
read_center_normal_from_tag( core::Vector & center, core::Vector & normal, utility::tag::TagCOP tag );

/// @brief Read in a user provided center/normal pair from the commandline, safetly
/// @details Read the membrane setup center & normal options form the command line
/// from mp:setup:center and mp:setup:normal. Intended to reduce IO code duplication
/// in membrane framework movers.
void
read_center_normal_from_cmd( core::Vector & center, core::Vector & normal );

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_util_hh

