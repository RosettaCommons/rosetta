// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/util.cc
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
///
/// NOTE: All of these methods require a RosettaMP framework pose or eventually
/// may require this. Use pose.conformation().is_membrane() for safety checks!
///
/// Last Modified: 7/9/15
/// @author Rebecca faye Alford (rfalford12@gmail.com)
/// @author JKLeman (julia.koehler1982@gmail.com)


// Unit Headers
#include <protocols/membrane/util.hh>

// Project Headers
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/Embedding.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/SpanningTopology.hh>

#include <core/conformation/membrane/Exceptions.hh>

// Package Headers
#include <protocols/moves/DsspMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>

#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/tag/Tag.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>

#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

#include <basic/Tracer.hh>

// C++ Headers
#include <algorithm>
#include <string>
#include <cstdlib>
#include <cmath>

static THREAD_LOCAL basic::Tracer TR( "protocols.membrane.util" );

using basic::Error;
using basic::Warning;

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
) {

	using namespace core::conformation::membrane;
	using namespace core::scoring;

	// Check that pose is actually a membrane pose
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot calculate membrane RMSD on a non membrane pose!" );
	}

	// Pick transmembrane spanning regions
	SpanningTopologyOP topology( pose.conformation().membrane_info()->spanning_topology() );
	ObjexxFCL::FArray1D_bool tm_regions ( pose.total_residue(), false );
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( topology->in_span(i) ) {
			tm_regions(i)=true;
		}
	}

	return rmsd_no_super_subset( native_pose, pose, tm_regions, is_protein_backbone );

}

/// @brief Compute all-atom RMSD between TM regions - don't superimpose
/// @details Calculate the rmsd between all atoms in the pose in the
/// transmembrane regions, as defined by the spanning topology object.
/// Do not superimpose the poses. Takes a native pose & current pose
core::Real
mem_all_atom_rmsd_no_super(
	core::pose::Pose & native_pose,
	core::pose::Pose & pose
) {

	using namespace core::conformation::membrane;
	using namespace core::scoring;

	// Check that pose is actually a membrane pose
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot calculate membrane RMSD on a non membrane pose!" );
	}

	// Pick transmembrane spanning regions
	core::conformation::membrane::SpanningTopologyOP topology( pose.conformation().membrane_info()->spanning_topology() );
	ObjexxFCL::FArray1D_bool tm_regions ( pose.total_residue(), false );
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( topology->in_span(i) ) {
			tm_regions(i)=true;
		}
	}

	return rmsd_no_super_subset( native_pose, pose, tm_regions, is_heavyatom );

}

/// @brief Compute backbone RMSD between TM regions - do superimpose
/// @details Calculate the rmsd between backbone atoms (N, CB, CA, O)
/// in the transmembrane regions, as defined by the spanning topology
/// object Superimpose the poses. Takes a native pose and current pose
core::Real
mem_bb_rmsd_with_super(
	core::pose::Pose & native_pose,
	core::pose::Pose & pose
) {

	using namespace core::conformation::membrane;
	using namespace core::scoring;

	// Check that pose is actually a membrane pose
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot calculate membrane RMSD on a non membrane pose!" );
	}

	// Pick transmembrane spanning regions
	SpanningTopologyOP topology( pose.conformation().membrane_info()->spanning_topology() );
	ObjexxFCL::FArray1D_bool tm_regions ( pose.total_residue(), false );
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( topology->in_span(i) ) {
			tm_regions(i)=true;
		}
	}

	return rmsd_with_super_subset( native_pose, pose, tm_regions, is_protein_backbone );

}

/// @brief Compute all-atom RMSD between TM regions - do superimpose
/// @details Calculate the rmsd between all atoms in the pose in the
/// transmembrane regions, as defined by the spanning topology object.
/// Superimpose the poses. Takes a native pose & current pose
core::Real
mem_all_atom_rmsd_with_super(
	core::pose::Pose & native_pose,
	core::pose::Pose & pose
)  {

	using namespace core::conformation::membrane;
	using namespace core::scoring;

	// Check that pose is actually a membrane pose
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot calculate membrane RMSD on a non membrane pose!" );
	}

	// Pick transmembrane spanning regions
	SpanningTopologyOP topology( pose.conformation().membrane_info()->spanning_topology() );
	ObjexxFCL::FArray1D_bool tm_regions ( pose.total_residue(), false );
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( topology->in_span(i) ) {
			tm_regions(i)=true;
		}
	}

	return rmsd_with_super_subset( native_pose, pose, tm_regions, is_heavyatom );

}

//////////////////////////////////////////////////////////////////
// Methods for calculating tilt of helices relative to membrane //
//////////////////////////////////////////////////////////////////

/// @brief Calculate tilt of a TM span relative to the membrane normal
/// @details Given a transmembrane span #, calculate the angle between the
/// axis through the helix and the membrane normal. Works for relatively
/// straight helices but less accurate for kinks. Takes a pose & span number
core::Real
calc_helix_tilt_angle(
	core::pose::Pose & pose,
	core::Size span_no
) {

	using namespace core;
	using namespace core::conformation::membrane;

	// Get membrane normal from membrane info
	Vector normal( pose.conformation().membrane_info()->membrane_normal() );

	// Calculate an axis representing the "axis of the helix"
	Vector helix_axis( calc_helix_axis( pose, span_no ) );

	// Normalize both vectors (just incase
	normal.normalize();
	helix_axis.normalize();

	// Check that both vectors have a positive orientation
	if ( normal.z() < 0 ) normal.negate();
	if ( helix_axis.z() < 0 ) helix_axis.negate();

	// Calculate the angle between the normal & helix
	core::Real tilt_angle( angle_of( normal, helix_axis ) );
	return numeric::conversions::degrees( tilt_angle );
}

/// @brief Determine the axis used to define a single TM Helix
/// @details Using the COM of the helix start & end position, calculate a helix
/// describing its geometry relative to the memrbane normal. Takes a pose &
/// span number. Not a good approx for helices with kinks.
/// TODO: CODE INSIDE SHOULD BE REPLACED WITH SPAN EMBEDDING CALCULATIONS!
///       THAT'S WHAT THEY ARE THERE FOR! (JKLeman)
core::Vector
calc_helix_axis(
	core::pose::Pose & pose,
	core::Size span_no
) {

	using namespace core;
	using namespace core::conformation::membrane;

	// Get the span from the pose
	SpanOP helix_span( pose.conformation().membrane_info()->spanning_topology()->span( span_no ) );

	// Check the core::Size of the span is sufficient for the this calculation
	if ( helix_span->end() - helix_span->start() < 6 ) {
		utility_exit_with_message( "Transmembrane span is too small to calculate a helix axis using the center of masses - less than six residues" );
	}

	// Check the core::Size of the span relative to the core::Size of the protein is appropriate for
	// a COM centered at the span start and end
	bool use_centered( true );
	if ( ( helix_span->start() < 2 ) ||
			( pose.total_residue() - helix_span->end() < 2 ) ) {
		use_centered = false;
	}

	// Grab CA coordinates from the pose based on specified method
	core::Vector start_0, start_1, start_2;
	core::Vector end_0, end_1, end_2;
	if ( use_centered ) {

		start_0 = pose.residue( helix_span->start() - 1 ).atom( "CA" ).xyz();
		start_1 = pose.residue( helix_span->start() ).atom( "CA" ).xyz();
		start_2 = pose.residue( helix_span->start() + 1 ).atom( "CA" ).xyz();

		end_0 = pose.residue( helix_span->end() - 1 ).atom( "CA" ).xyz();
		end_1 = pose.residue( helix_span->end() ).atom( "CA" ).xyz();
		end_2 = pose.residue( helix_span->end() + 1 ).atom( "CA" ).xyz();

	} else {

		start_0 = pose.residue( helix_span->start() ).atom( "CA" ).xyz();
		start_1 = pose.residue( helix_span->start() + 1 ).atom( "CA" ).xyz();
		start_2 = pose.residue( helix_span->start() + 2 ).atom( "CA" ).xyz();

		end_0 = pose.residue( helix_span->end() ).atom( "CA" ).xyz();
		end_1 = pose.residue( helix_span->end() -1 ).atom( "CA" ).xyz();
		end_2 = pose.residue( helix_span->end() -2 ).atom( "CA" ).xyz();

	}

	// Calculate the vector as the difference in start & end com
	core::Vector helix_axis = com( end_0, end_1, end_2 ) - com( start_0, start_1, start_2 );
	return helix_axis;

}

/// @brief Calculate center of mass between 3 xyz coords
/// @details Given three xyz vectors, calculate the center of mass
/// and return a vector. Helper method to calc_helix axis.
core::Vector
com( core::Vector a, core::Vector b, core::Vector c ) {

	core::Real x( (a.x() + b.x() + c.x()) / 3 );
	core::Real y( (a.y() + b.y() + c.y()) / 3 );
	core::Real z( (a.z() + b.z() + c.z()) / 3 );
	return core::Vector( x, y, z);

}

/// @brief Calculate the RMSD between a helix tilt angle & reference
/// @details Given a reference angle and measured angle, calculate the
/// root mean square deviation between the two single values. Takes
/// the measured tilt angle and reference angle (typically from experiment)
core::Real
calc_angle_rmsd( core::Real measured_angle, core::Real ref_angle ) {

	core::Real abs_diff( std::abs( measured_angle - ref_angle ) );
	core::Real rms = std::sqrt( ( std::pow( abs_diff, 2 ) / 2 ) );
	return rms;
}

/// @brief Calculate tilt angle and distance from membrane center
/// @details Computes the tilt angle and distance from membrane center for the
///   protein embedding center and normal
utility::vector1< core::Real > pose_tilt_angle_and_center_distance( core::pose::Pose & pose ) {

	using namespace protocols::membrane::geometry;
	utility::vector1< core::Real > angle_and_distance;

	// membrane center
	core::Vector mem_center = pose.conformation().membrane_info()->membrane_center();
	core::Vector mem_normal = pose.conformation().membrane_info()->membrane_normal();

	// compute structure-based embedding
	EmbeddingDefOP emb( compute_structure_based_embedding( pose ) );

	// compute tilt angle between embedding normal and membrane normal
	core::Real tilt_angle( numeric::conversions::degrees( angle_of( mem_normal, emb->normal() ) ) );
	TR << "tilt angle: " << tilt_angle << std::endl;

	// compute distance from membrane center plane
	core::Real dist_from_center = dot( emb->center() - mem_center, mem_normal );
	TR << "dist from center plane: " << dist_from_center << std::endl;

	angle_and_distance.push_back( tilt_angle );
	angle_and_distance.push_back( dist_from_center );

	return angle_and_distance;

} // compute tilt angle and distance from membrane center

////////////////////////////////////////////////////////////////
// Safety checks & convenience methods for membrane foldtrees //
////////////////////////////////////////////////////////////////

/// @brief Determine whether the membrane is modeled as fixed
/// @details Based on the setup of the foldtree, determined whether
/// the membrane is currently fixed, meaning it is setup as the
/// root in the FoldTree and has no upstream children. Takes a pose.
bool
is_membrane_fixed( core::pose::Pose & pose ) {

	using namespace core::conformation::membrane;
	using namespace core::kinematics;

	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message("Pose is not a membrane pose. Quitting.");
	}

	// Get membrane res, jump & upstream residue
	core::Size membrane_rsd( pose.conformation().membrane_info()->membrane_rsd_num() );
	core::Size membrane_jump( pose.conformation().membrane_info()->membrane_jump() );
	core::Size upstream_res( pose.conformation().fold_tree().upstream_jump_residue( membrane_jump ) );

	if ( upstream_res == membrane_rsd &&
			pose.conformation().fold_tree().is_root( int( upstream_res ) )
			) {
		return true;
	}

	return false;
}

/// @brief Determine whether membrane can move on its own
/// @details Based on the setup of the FoldTree, determine whether
/// the membrane is moveable, but when moved, won't cause anything in the
/// protein to move (i.e. independently moveable). Takes a pose.
bool
is_membrane_moveable_by_itself( core::pose::Pose & pose ) {

	using namespace core::kinematics;

	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message("Pose is not a membrane pose. Quitting.");
	}

	// If fixed, return false
	if ( is_membrane_fixed( pose ) ) return false;

	// Grab the current foldtree from the conformation
	FoldTree const & current_ft( pose.conformation().fold_tree() );

	// Grab membrane info from the pose
	core::Size membrane_rsd( pose.conformation().membrane_info()->membrane_rsd_num() );
	core::Size membrane_jump( pose.conformation().membrane_info()->membrane_jump() );

	if ( current_ft.num_jump() > 1 ) {

		// Iterate through the edge list and check that only one jump
		// connects to the membrane residue
		for ( core::kinematics::FoldTree::const_iterator it = current_ft.begin(), it_end = current_ft.end(); it != it_end; ++it ) {

			// If a an edge is a jump edge that is not the membrane jump
			// but has a start or end point that is the membrane rsd,
			// the memrbane rsd is not 'independently moveable'
			if ( ( it->label() > 0 ) && (
					( it->label() != int(membrane_jump) ) &&
					( it->start() == int(membrane_rsd) ||
					it->stop() == int(membrane_rsd) ) ) ) {
				return false;
			}
		}
	}

	return true;
}

/// @brief Set membrane residue to root of foldtree
/// @details Naively sets the root of the foldtree to be the membrane
/// residue. Should perform checks before doing this!
void reorder_membrane_foldtree( core::pose::Pose & pose ) {

	// get foldtree from pose
	core::kinematics::FoldTree foldtree = pose.fold_tree();

	// reorder foldtree
	foldtree.reorder( pose.conformation().membrane_info()->membrane_rsd_num() );

	// set foldtree in pose
	pose.fold_tree( foldtree );
} // reorder membrane foldtree

//////////////////////////////////////////////////////////////

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
///  iJ = interface jump, will be returned from the function
///
core::Size create_membrane_docking_foldtree_from_partners( core::pose::Pose & pose, std::string const partners ) {

	using namespace utility;
	using namespace core;
	using namespace core::kinematics;

	// check that partners isn't empty
	if ( partners.size() == 0 || partners == "_" ) {
		utility_exit_with_message( "Membrane docking partners are undefined. Quitting..." );
	}

	// split partner string (AB, CDE)
	utility::vector1< std::string > partner( utility::string_split( partners, '_' ) );

	// initialize partners with chains (will be 1,2 / 3,4,5)
	// initialize anchor points within these chains (will be 19,234 / 287,354,528)
	// (the anchor points are the chain TM COMs)
	utility::vector1< core::Size > chains1;
	utility::vector1< core::Size > chains2;
	utility::vector1< core::Size > anchors1;
	utility::vector1< core::Size > anchors2;
	utility::vector1< core::Size > cutpoints1;
	utility::vector1< core::Size > cutpoints2;

	// go through first partner chainIDs, convert into chain numbers, add to vector
	// also get anchor points and cutpoints for these chains
	for ( core::Size i = 1; i <= partner[ 1 ].size(); ++i ) {

		// get chain, add to chains vector and get anchor point
		core::Size chain = get_chain_id_from_chain( partner[ 1 ][ i-1 ], pose );
		chains1.push_back( chain );
		anchors1.push_back( rsd_closest_to_chain_tm_com( pose, chain ) );
		cutpoints1.push_back( chain_end_res( pose, chain ) );
	}

	// go through second partner chainIDs, convert into chain numbers, add to vector
	// also get anchor points and cutpoints for these chains
	for ( core::Size i = 1; i <= partner[ 2 ].size(); ++i ) {

		// get chain, add to chains vector and get anchor point
		core::Size chain = get_chain_id_from_chain( partner[ 2 ][ i-1 ], pose );
		chains2.push_back( chain );
		anchors2.push_back( rsd_closest_to_chain_tm_com( pose, chain ) );
		cutpoints2.push_back( chain_end_res( pose, chain ) );
	}

	// create simple foldtree
	FoldTree ft = FoldTree();
	ft.simple_tree( pose.total_residue() );

	// get membrane residue
	core::Size memrsd = pose.conformation().membrane_info()->membrane_rsd_num();

	TR << "partners " << partners << std::endl;
	TR << "mem rsd " << memrsd << std::endl;
	TR << "anchors[1] " << anchors1[1] << std::endl;


	// anchor MEM on the first chain of the first partner with the cutpoint
	// right before the MEM residues
	ft.new_jump( memrsd, anchors1[ 1 ], memrsd - 1 );

	// create jumps between the chains in partner1
	for ( core::Size i = 2; i <= anchors1.size(); ++i ) {
		ft.new_jump( anchors1[ 1 ], anchors1[ i ], cutpoints1[ i-1 ] );
	}

	// create jumps between the chains in partner1
	for ( core::Size i = 2; i <= anchors2.size(); ++i ) {
		ft.new_jump( anchors2[ 1 ], anchors2[ i ], cutpoints2[ i-1 ] );
	}

	// create interface jump between the partners by connecting their 1st chains
	// cutpoint is cutpoint of last chain in partner 1
	int interface_jump = ft.new_jump( anchors1[ 1 ], anchors2[ 1 ], cutpoints1[ cutpoints1.size() ] );

	// reorder and set the foldtree in the pose to the newly created one
	ft.reorder( memrsd );
	ft.show( TR );
	pose.fold_tree( ft );

	// set the membrane jump in MembraneInfo
	pose.conformation().membrane_info()->set_membrane_jump( static_cast< core::SSize >( 1 ) );

	return static_cast< core::Size >( interface_jump );

} // create_membrane_docking_foldtree_from_partners

////////////////////////////////////////////////////////////////////////////////

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
core::Size create_membrane_foldtree_anchor_com( core::pose::Pose & pose ) {

	using namespace core::kinematics;
	using namespace core::pose;

	// if pose not membrane pose, cry
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Can't create a membrane foldtree on a non-membrane pose. Quitting..." );
	}

	// get anchor points for jumps: get all chainids
	utility::vector1< int > chains = get_chains( pose );

	// initialize empty vector for anchor points (i.e. jump rsd positions)
	utility::vector1< core::Size > anchors;

	// get residues closest to COMs for all chains which will be new jump anchor residues
	for ( core::Size i = 1; i < chains.size(); ++i ) {
		core::Size anchor = rsd_closest_to_chain_com( pose, chains[ i ] );
		anchors.push_back( anchor );
	}

	// create foldtree
	create_membrane_foldtree_from_anchors( pose, anchors );

	return anchors[1];

} // create membrane foldtree anchor center-of-mass

/////////////////////////////////////////

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
core::Size create_membrane_foldtree_anchor_tmcom( core::pose::Pose & pose ) {

	using namespace core::kinematics;

	// if pose not membrane pose, cry
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Can't create a membrane foldtree on a non-membrane pose. Quitting..." );
	}

	// get the chain TM COM anchor points
	utility::vector1< core::Size > anchors( get_anchor_points_for_tmcom( pose ) );

	// create foldtree
	create_membrane_foldtree_from_anchors( pose, anchors );

	return anchors[1];

} // create membrane foldtree anchor tm center-of-mass

/////////////////////////////////////////

/// @brief Create membrane foldtree from scratch
/// @details The foldtree is setup such that the membrane is at the root and
///   anchored at the residue closest to the pose TM COM with jumps from there
///   to each chain TRANSMEMBRANE COM; requires the membrane to be present
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
core::Size create_membrane_foldtree_anchor_pose_tmcom( core::pose::Pose & pose ) {

	using namespace core::kinematics;
	using namespace core::pose;

	// if pose not membrane pose, cry
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Can't create a membrane foldtree on a non-membrane pose. Quitting..." );
	}

	// get anchor points for jumps: get all chainids
	utility::vector1< int > chains = get_chains( pose );

	// initialize empty vector for anchor points (i.e. jump rsd positions)
	utility::vector1< core::Size > anchors;

	// get rsd nearest pose TM COM
	core::Size root_anchor( rsd_closest_to_pose_tm_com( pose ) );

	// get chain for pose TM COM
	core::Size root_anchor_chain( static_cast< core::Size >( pose.chain(root_anchor) ) );

	// push back root anchor into anchors vector
	anchors.push_back( root_anchor );

	// get residues closest to COMs for all chains which will be new jump anchor residues
	for ( core::Size i = 1; i < chains.size(); ++i ) {

		// root anchor chain is already in vector at the first position
		if ( i != root_anchor_chain ) {
			core::Size anchor = rsd_closest_to_chain_tm_com( pose, chains[ i ] );
			anchors.push_back( anchor );
		}
	}

	// debug
	// for ( core::Size i = 1; i <= anchors.size(); ++i ) {
	//  TR << "anchor " << i << ": " << anchors[i] << std::endl;
	// }

	// create foldtree
	create_membrane_foldtree_from_anchors( pose, anchors );

	return anchors[1];

} // create_membrane_foldtree_anchor_pose_tmcom

/////////////////////////////////////////

/// @brief Create foldtree for multiple docking partners
utility::vector1< core::Size > create_membrane_multi_partner_foldtree_anchor_tmcom( core::pose::Pose & pose, std::string partner ) {

	// get chains from pose
	utility::vector1< int > chainids( get_chains( pose ) );
	utility::vector1< core::Size > chain_ends( chain_end_res( pose ) );

	// get chain TM COMs as anchor points
	utility::vector1< core::Vector > all_anchors;
	utility::vector1< core::Size > anchors ( get_anchor_points_for_tmcom( pose ) );

	// print anchors
	for ( core::Size i = 1; i <= anchors.size(); ++i ) {
		TR << "anchor " << i << " " << anchors[ i ] << std::endl;
	}

	// add jump from MEM to first chain anchor
	core::Size memrsd = pose.conformation().membrane_info()->membrane_rsd_num();
	core::Size cutpoint = chain_ends[ 1 ];
	core::Vector jump1( anchors[ 1 ], memrsd, cutpoint );
	all_anchors.push_back( jump1 );
	core::Size jumpnumber( 1 );
	utility::vector1< core::Size > jumps;

	// split partners by underscore
	utility::vector1< std::string > partners( utility::string_split( partner, '_' ) );

	// go through partners (AB / CD)
	for ( core::Size p = 1; p <= partners.size(); ++p ) {

		// get anchor for first chain
		core::Size chainid;
		core::Size first_chain_anchor( 0 );

		// convert chain from partners (f.ex. B) into chainID (f.ex. 2)
		chainid = get_chain_id_from_chain( partners[ p ][ 0 ], pose );

		// go through chainids (except MEM) to find chain of interest  and first anchor point
		for ( core::Size i = 1; i <= chainids.size()-1; ++i ) {
			if ( chainid == static_cast< core::Size >( chainids[ i ] ) ) {
				first_chain_anchor = anchors[ i ];
			}
		}

		// add jump from first chain to this chain
		if ( first_chain_anchor != anchors[ 1 ] ) {
			++jumpnumber;
			jumps.push_back( jumpnumber );
			cutpoint = chain_ends[ chainid ];
			core::Vector jump( anchors[ 1 ], first_chain_anchor, cutpoint );
			all_anchors.push_back( jump );
		}

		// if multiple chains: anchor first chain on first chain, then the next
		// on the first chain in the partner
		if ( partners[ p ].size() > 1 ) {

			core::Size chain_anchor( 0 );

			// iterate through rest of chains in the partner (1,2,3), string indexes with 0
			for ( core::Size c = 1; c < partners[ p ].size(); ++c ) {

				// go through chainids (except MEM) to find chain of interest and anchor point
				for ( core::Size i = 1; i <= chainids.size()-1; ++i ) {

					// convert chain from partners (f.ex. B) into chainID (f.ex. 2)
					chainid = get_chain_id_from_chain( partners[ p ][ c ], pose );

					// get anchor point of current partner
					if ( chainid == static_cast< core::Size >( chainids[ i ] ) ) {
						chain_anchor = anchors[ chainid ];
					}

				} // iterate over chainids

				// add jump from first chain anchor in the partner to this anchor point
				cutpoint = chain_ends[ chainid ];
				core::Vector jump( first_chain_anchor, chain_anchor, cutpoint );
				all_anchors.push_back( jump );

			} // iterate through rest of chains in the partner

		} // multiple chains in partner

	} // go through partners

	for ( core::Size i = 1; i <= all_anchors.size(); ++i ) {
		TR << all_anchors[ i ].to_string() << std::endl;
	}

	// create foldtree from anchors
	create_specific_membrane_foldtree( pose, all_anchors );

	return jumps;

} // create membrane multi partner foldtree anchor tmcom

/////////////////////////////////////////

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
void create_membrane_foldtree_from_anchors(
	core::pose::Pose & pose,
	utility::vector1< core::Size > anchors
) {

	using namespace core::kinematics;
	using namespace core::pose;

	// if pose not membrane pose, cry
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Can't create a membrane foldtree on a non-membrane pose. Quitting..." );
	}

	// [1]
	// get all chain poses to get number of residues in each chain
	// this is to define the cutpoints later on
	utility::vector1< PoseOP > chain_poses = pose.split_by_chain();

	// go through chains to get number or residues in each chain
	// counter needs to be smaller than the number of chains because the
	// membrane residue has its own chain at the end
	utility::vector1< core::Size > chain_nres;
	core::Size rsd_counter( 0 );
	for ( core::Size i = 1; i < chain_poses.size(); ++i ) {
		rsd_counter += chain_poses[ i ]->total_residue();
		chain_nres.push_back( rsd_counter );
	}

	// [2]
	// get membrane residue
	core::Size memrsd = pose.conformation().membrane_info()->membrane_rsd_num();

	// [3]
	// start with simple tree
	FoldTree ft = FoldTree();
	ft.simple_tree( pose.total_residue() );

	// add jumps with cutpoints, returns jump number of new jump:
	// add first jump from membrane residue to first chain
	ft.new_jump( memrsd, anchors[ 1 ], memrsd - 1 );

	// [4]
	// for each chain starting with the second one, add another jump with
	// the cutpoint after the chain
	for ( core::Size i = 2; i <= chain_poses.size()-1; ++i ) {
		TR << "adding jump from " << anchors[1] << " to " << anchors[ i ] << " with cutpoint " << chain_nres[ i-1 ] << std::endl;
		ft.new_jump( anchors[ 1 ], anchors[ i ], chain_nres[ i-1 ] );
	}

	// set the foldtree in the pose to the newly created one
	ft.reorder( memrsd );
	ft.show( TR );
	pose.fold_tree( ft );

	// set the membrane jump in MembraneInfo
	pose.conformation().membrane_info()->set_membrane_jump( static_cast< core::SSize >( 1 ) );

} // create membrane foldtree from anchors

////////////////////////////////////////////////////////////////////////////////

/// @brief Helper function to create a specific membrane foldtree
/// @details I am hijacking xyzVectors to hold the jumps that need to be
///   created in the foldtree: xyz = ( rsd1, rsd2, cutpoint )
///   THE JUMP NUMBERS WILL BE CONSECUTIVE, ACCORDING TO THE VECTOR1
core::Size create_specific_membrane_foldtree( core::pose::Pose & pose, utility::vector1< core::Vector > anchors ) {

	using namespace core::kinematics;

	// if pose not membrane pose, cry
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Can't create a membrane foldtree on a non-membrane pose. Quitting..." );
	}

	// create a simple tree
	FoldTree ft = FoldTree();
	ft.simple_tree( pose.total_residue() );

	// membrane residue and jump
	core::Size mem_rsd = pose.conformation().membrane_info()->membrane_rsd_num();
	core::Size mem_jump = static_cast< core::Size > ( pose.conformation().membrane_info()->membrane_jump() );

	// go through anchor points
	for ( core::Size i = 1; i <= anchors.size(); ++i ) {

		TR << "adding jump from " << anchors[i].x() << " to " << anchors[i].y()
			<< " with cutpoint " << anchors[i].z() << std::endl;

		// add jumps with cutpoints, returns jump number of new jump:
		ft.new_jump( anchors[ i ].x(), anchors[ i ].y(), anchors[ i ].z() );

		// find membrane jump
		if ( anchors[ i ].x() == mem_rsd || anchors[ i ].y() == mem_rsd ) {
			mem_jump = i;
		}
	}

	// set the foldtree in the pose to the newly created one
	ft.reorder( mem_rsd );
	ft.show( TR );
	pose.fold_tree( ft );

	// set the membrane jump in MembraneInfo
	pose.conformation().membrane_info()->set_membrane_jump( static_cast< core::SSize >( mem_jump ) );

	// returns main anchor for MEM on first chain
	return anchors[1].x();

} // create specific membrane foldtree

////////////////////////////////////////////////////////////////////////////////

/// @brief Helper function to create membrane foldtrees
/// @details Returns the residues closest to the COMs for each chain
utility::vector1< core::Size > get_anchor_points_for_tmcom( core::pose::Pose & pose ) {

	// get anchor points for jumps: get all chainids
	utility::vector1< int > chains = get_chains( pose );

	// initialize empty vector for anchor points (i.e. jump rsd positions)
	utility::vector1< core::Size > anchors;

	// get residues closest to COMs for all chains which will be new jump anchor residues
	// needs to < chains.size() because the MEM is an additional chain
	for ( core::Size i = 1; i < chains.size(); ++i ) {
		TR << "chain " << i << " " << chains[ i ] << std::endl;
		core::Size anchor = protocols::membrane::rsd_closest_to_chain_tm_com( pose, chains[ i ] );
		anchors.push_back( anchor );
	}

	return anchors;
} //  get anchor points for TM COM

/////////////////////////////////////////

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
core::Size setup_foldtree_pose_com( core::pose::Pose & pose ) {

	using namespace core::kinematics;
	using namespace core::pose;

	// if pose is membrane pose, cry
	if ( pose.conformation().is_membrane() ) {
		utility_exit_with_message( "There are better foldtrees for membrane proteins than this: check out protocols/membrane/util. Quitting." );
	}

	// get anchor points for jumps: get all chainids
	utility::vector1< int > chains = get_chains( pose );

	// initialize empty vector for anchor points (i.e. jump rsd positions)
	utility::vector1< core::Size > anchors;

	// get rsd nearest pose COM
	core::Size root_anchor( static_cast< core::Size >( residue_center_of_mass( pose, 1, pose.total_residue() ) ) );

	// get chain for pose COM
	core::Size root_anchor_chain( static_cast< core::Size >( pose.chain(root_anchor) ) );

	// push back root anchor into anchors vector
	anchors.push_back( root_anchor );

	// get residues closest to COMs for all chains which will be new jump anchor residues
	for ( core::Size i = 1; i < chains.size(); ++i ) {

		// root anchor chain is already in vector at the first position
		if ( i != root_anchor_chain ) {
			core::Size anchor = rsd_closest_to_chain_com( pose, chains[ i ] );
			anchors.push_back( anchor );
		}
	}

	// setup foldtree
	setup_foldtree_from_anchors( pose, anchors );

	return anchors[1];

} // setup foldtree chain COM

/////////////////////////////////////////

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
void setup_foldtree_from_anchors( core::pose::Pose & pose, utility::vector1< core::Size > anchors ) {

	using namespace core::kinematics;
	using namespace core::pose;

	// if pose is membrane pose, cry
	if ( pose.conformation().is_membrane() ) {
		utility_exit_with_message( "There are better foldtrees for membrane proteins than this: check out protocols/membrane/util. Quitting." );
	}

	// get all chain poses to get number of residues in each chain
	// this is to define the cutpoints later on
	utility::vector1< PoseOP > chain_poses = pose.split_by_chain();

	// go through chains to get number or residues in each chain
	utility::vector1< core::Size > chain_nres;
	core::Size rsd_counter( 0 );
	for ( core::Size i = 1; i <= chain_poses.size(); ++i ) {
		rsd_counter += chain_poses[ i ]->total_residue();
		chain_nres.push_back( rsd_counter );
	}

	// start with simple tree
	FoldTree ft = FoldTree();
	ft.simple_tree( pose.total_residue() );

	// for each chain starting with the second one, add another jump with
	// the cutpoint after the chain
	for ( core::Size i = 2; i <= chain_poses.size()-1; ++i ) {
		ft.new_jump( anchors[ 1 ], anchors[ i ], chain_nres[ i-1 ] );
	}

	// reorder the foldtree such that the first anchor point is at the root
	ft.reorder( anchors[ 1 ] );

	// set the foldtree in the pose to the newly created one
	ft.show( TR );
	pose.fold_tree( ft );

} // setup foldtree from anchors

///////////////////////////////////////////////////////////
// Utilities for accessing dssp, z coords and chain info //
///////////////////////////////////////////////////////////

/// @brief Grab the z-coordinates and chainIDs from the entire pose
/// @details From the pose, grab all of the z_coords of CA atoms and
/// chain IDs, currently used for spanning topology construction.
/// Returns a std::pair of two vectors: the first a vector1 of z
/// coordinates and the second a vector1 of chainIDs for CA atoms
std::pair< utility::vector1< core::Real >, utility::vector1< core::Real > >
get_chain_and_z( core::pose::Pose const & pose ) {

	TR.Debug << "get_pose_info" << std::endl;
	using namespace core::pose;

	// initialize variables
	utility::vector1< core::Real > z_coord;
	utility::vector1< core::Size > chain_info;

	// loop over residues, get chain info and z_coord
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {

		// get info
		z_coord.push_back( static_cast< core::Real >(pose.residue(i).atom(2).xyz().z()) );
		chain_info.push_back( pose.chain(i) );
	}

	// put the data in a pair
	std::pair< utility::vector1< core::Real >, utility::vector1< core::Size > > pose_info( z_coord, chain_info );

	return pose_info;

} // get chain and z from pose

/// @brief  Get dssp defined secondary structure from the pose
/// @details Given a pose, grab a vector of characters describing the secondary
/// structure at each residue position in the pose, defined by DSSP
utility::vector1< char > get_secstruct( core::pose::Pose & pose ) {

	TR.Debug << "get_secstruct" << std::endl;
	using namespace core::pose;

	// set secondary structure in pose with DSSP
	protocols::moves::DsspMover dssp = protocols::moves::DsspMover();
	dssp.apply( pose );

	// initialize variable
	utility::vector1< char > secstruct;

	// loop over residues, get chain info and z_coord
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {

		// get info
		secstruct.push_back( pose.conformation().secstruct( i ) );
	}

	return secstruct;

}// get_secstruct


///////////////////////////////////////////////////////////////////
// Methods for calculating the protein embedding in the membrane //
///////////////////////////////////////////////////////////////////

/// @brief Compute Membrane Center/Normal from Membrane Spanning
/// topology
void compute_structure_based_embedding(
	core::pose::Pose const & pose,
	core::conformation::membrane::SpanningTopology const & topology,
	core::Vector & center,
	core::Vector & normal
) {

	using namespace protocols::membrane::geometry;

	// create EmbeddingDef to return
	EmbeddingDefOP embed = compute_structure_based_embedding( pose, topology );

	// set new center and normal
	center = embed->center();
	normal = embed->normal();

}// compute structure-based position


/// @brief Compute Membrane Center/Normal from Membrane Spanning
/// topology, uses topology from MembraneInfo
void compute_structure_based_embedding(
	core::pose::Pose const & pose,
	core::Vector & center,
	core::Vector & normal
) {

	using namespace core::conformation::membrane;

	// get topology from MembraneInfo
	SpanningTopology topo( *pose.conformation().membrane_info()->spanning_topology() );

	// create EmbeddingDef to return
	compute_structure_based_embedding( pose, topo, center, normal );

}// compute structure-based position

/// @brief Compute Membrane Center/Normal from Membrane Spanning
/// topology
protocols::membrane::geometry::EmbeddingDefOP
compute_structure_based_embedding(
	core::pose::Pose const & pose,
	core::conformation::membrane::SpanningTopology const & topo
) {

	using namespace protocols::membrane::geometry;
	using namespace core::conformation::membrane;

	if ( topo.nspans() == 0 ) {
		utility_exit_with_message("The SpanningTopology object in MembraneInfo is empty!" );
	}

	// create Embedding object
	Embedding embeddings = Embedding( topo, pose );

	// return total embedding
	return embeddings.total_embed();

}// compute structure based membrane position


/// @brief Compute Membrane Center/Normal from Membrane Spanning
/// topology, uses topology from MembraneInfo
protocols::membrane::geometry::EmbeddingDefOP
compute_structure_based_embedding( core::pose::Pose const & pose ){

	using namespace core::conformation::membrane;

	// get topology from MembraneInfo
	SpanningTopology topo( *pose.conformation().membrane_info()->spanning_topology() );

	return compute_structure_based_embedding( pose, topo );

}// compute structure based membrane position

/// @brief Compute embedding by chain
/// @details The embeddings can be computed either from pose and topology or they
///   can be optimized differently; the function correlates each EmbeddingDef
///   object in embeddings with a span object in the pose's topology
protocols::membrane::geometry::EmbeddingOP
compute_embeddings_by_chain( core::pose::Pose const & pose ) {

	using namespace core::conformation::membrane;
	using namespace protocols::membrane::geometry;

	// get topology from pose
	SpanningTopologyOP topo = pose.conformation().membrane_info()->spanning_topology();

	// split topology by chain
	utility::vector1< SpanningTopologyOP > chain_topos( split_topology_by_chain_noshift( pose, topo ) );

	// initialize vector of EmbeddingDefs
	EmbeddingOP embeddings( new Embedding() );

	// go through each chain
	for ( core::Size i = 1; i <= chain_topos.size(); ++i ) {

		// compute structure-based embedding for this spanning topology
		EmbeddingDefOP chain_emb = compute_structure_based_embedding( pose, *chain_topos[ i ] );

		// push back in vector
		embeddings->add_span_embedding( chain_emb );
	}

	return embeddings;

} // compute embeddings by chain


/// @brief Average EmbeddingDefs as they are without vector inversion accounting for topology
/// @details Get average center and normal from a vector of EmbeddingDefs
protocols::membrane::geometry::EmbeddingDefOP
average_embeddings( utility::vector1< protocols::membrane::geometry::EmbeddingDefOP > const parts ) {

	using namespace protocols::membrane::geometry;

	// Initialize vars
	core::Vector center(0, 0, 0);
	core::Vector normal(0, 0, 0);

	// Compute resulting center and normal
	for ( core::Size i = 1 ; i <= parts.size(); ++i ) {
		center += parts[i]->center();
		normal += parts[i]->normal();
	}

	center /= parts.size();
	normal.normalize();

	// Create new embedding setup and return it
	EmbeddingDefOP embedding( new EmbeddingDef( center, normal ) );
	return embedding;

}// average embeddings

/// @brief Average EmbeddingDefs after first inverting some vectors accounting for topology
/// @details Get average center and normal from a vector of EmbeddingDefs
protocols::membrane::geometry::EmbeddingDefOP
average_antiparallel_embeddings(
	utility::vector1< protocols::membrane::geometry::EmbeddingDefOP > const parts
) {

	using namespace protocols::membrane::geometry;

	// Initialize vars
	core::Vector center(0, 0, 0);
	core::Vector normal(0, 0, 0);

	// embedding of first span
	core::Vector const center1 = parts[1]->center();
	core::Vector const normal1 = parts[1]->normal();

	// Compute resulting center and normal
	for ( core::Size i = 1 ; i <= parts.size(); ++i ) {

		TR << "center: " << parts[i]->center().to_string() << "normal: " << parts[i]->normal().to_string() << std::endl;

		// calculate new center
		center += parts[i]->center();

		// calculate points for angle calculation
		core::Vector p1 = center1 + normal1;
		core::Vector p  = center1 + parts[i]->normal();

		// calculate  angle between normals of first object and this one
		core::Real angle( numeric::angle_degrees( p1, center1, p ) );
		TR << "angle: " << angle << std::endl;

		// check if angle of normal is < 100 degrees to first normal
		// if yes, then add to normal, if no add inverted vector
		if ( angle > -100 && angle < 100 ) {
			normal += parts[i]->normal();
		} else {
			normal -= parts[i]->normal();
		}
	}

	center /= parts.size();
	normal.normalize();

	// Create new embedding setup and return it
	EmbeddingDefOP embedding( new EmbeddingDef( center, normal ) );
	return embedding;

}// average antiparallel embeddings

////////////////////////////////////////////////////////////////////////////////

/// @brief Update embedding of the partners after a move
/// @details Requires the jump number between the partners, the topology will
///    be taken from MembraneInfo and will be split accordingly; up and
///    down means upstream and downstream
void
update_partner_embeddings(
	core::pose::Pose const & pose,
	core::Size const jumpnum,
	protocols::membrane::geometry::EmbeddingDef & emb_up,
	protocols::membrane::geometry::EmbeddingDef & emb_down
) {

	using namespace protocols::membrane::geometry;
	using namespace core::conformation::membrane;

	// SpanningTopology objects
	SpanningTopologyOP topo = pose.conformation().membrane_info()->spanning_topology();
	SpanningTopologyOP topo_up( new SpanningTopology() ); // upstream partner
	SpanningTopologyOP topo_down( new SpanningTopology() ); // downstream partner

	// splitting topology by jump into upstream and downstream topology
	split_topology_by_jump_noshift( pose, jumpnum, topo, topo_up, topo_down );

	// compute embedding for partners (compute structure-based embedding with split topologies)
	EmbeddingDefOP emb1( compute_structure_based_embedding( pose, *topo_up ) );
	EmbeddingDefOP emb2( compute_structure_based_embedding( pose, *topo_down ) );

	// create new embedding objects to be able to dereference the pointer
	emb_up = EmbeddingDef( emb1->center(), emb1->normal() );
	emb_down = EmbeddingDef( emb2->center(), emb2->normal() );

} // update partner embeddings

////////////////////////////////////////////////////////////////////////////////

/// @brief Pose transmembrane center-of-mass
/// @details Gets the coordinates of the TM span center-of-mass
///   This only looks at the span start and end residues for
///   calculation of the TM span COM, this should be faster than the
///   real thing though
core::Vector pose_tm_com( core::pose::Pose const & pose ) {

	using namespace core::conformation::membrane;

	// get topology from MembraneInfo
	SpanningTopologyOP topo( pose.conformation().membrane_info()->spanning_topology() );

	// initialize vector
	core::Vector com( 0, 0, 0 );

	// go through topology and avg coords of start and end residues
	for ( core::Size i = 1; i <= topo->nspans(); ++i ) {

		// get CA coords
		core::Vector start_coord = pose.residue( topo->span(i)->start() ).xyz("CA");
		core::Vector end_coord = pose.residue( topo->span(i)->end() ).xyz("CA");

		// add to com vector
		com += start_coord;
		com += end_coord;
	}

	// divide by the number of points
	com /= ( 2 * topo->nspans() );

	return com;

} // pose TM COM

/////////////////////////////////////////

/// @brief Chain center-of-mass
/// @details Gets the coordinates of the chain center-of-mass
core::Vector
chain_com( core::pose::Pose const & pose, int chainid ) {

	using namespace core::pose;

	// split pose by chain
	PoseOP pose_chain( pose.split_by_chain( chainid ) );

	// compute COM of the chain
	core::Vector com( get_center_of_mass( *pose_chain ) );

	return com;

} // chain COM

////////////////////////////////////////////////////////////////////////////////

/// @brief Chain center-of-mass of TM regions
/// @details Gets the coordinates of the chain center-of-mass but only the TM regions
///   BE AWARE THAT THE LAST CHAIN FOR MEMBRANE PROTEINS IS THE MEMBRANE RESIDUE!!!
core::Vector chain_tm_com( core::pose::Pose const & pose, int chain ) {

	using namespace core::conformation::membrane;

	// check that the chain isn't the membrane residue
	if ( chain == pose.chain( pose.conformation().membrane_info()->membrane_rsd_num() ) ) {
		utility_exit_with_message("You are trying to compute the center-of-mass for the membrane residue as a chain. Choose a different one...");
	}

	// get topology and split it by chain without shifting the numbering in topology
	SpanningTopologyOP topo = pose.conformation().membrane_info()->spanning_topology();
	utility::vector1< SpanningTopologyOP > split_topo = split_topology_by_chain_noshift( pose, topo );

	// initializations
	utility::vector1< core::Size > splice_rsd;
	core::Size rsd_counter( 0 );

	// get vector of residues to take into account for COM calculation
	// go through all topologies
	for ( core::Size i = 1; i <= split_topo.size(); ++ i ) {

		// if we are looking at the topology of the chain of interest
		if ( i == static_cast< core::Size >( chain ) ) {

			// go through topology of the chain
			for ( core::Size j = 1; j <= split_topo[ i ]->nspans(); ++j ) {

				// set residue counter to first residue in the span
				rsd_counter = split_topo[ i ]->span( j )->start();

				// add to splice residues and increase until end of span
				while ( rsd_counter <= split_topo[ i ]->span( j )->end() ) {

					splice_rsd.push_back( rsd_counter );
					rsd_counter++;
				}
			}
		}
	}

	// if there are no splice residues because the topology is empty for that chain
	// i.e. for a soluble chain, get the COM of that chain without topology
	if ( splice_rsd.size() == 0 ) {
		core::pose::PoseOP sol_subpose = pose.split_by_chain( static_cast< core::Size >( chain ) );
		core::Vector sol_com( get_center_of_mass( *sol_subpose ) );
		return sol_com;
	}

	// debug
	// for ( core::Size i = 1; i <= splice_rsd.size(); ++i ) {
	//  TR << "splice res " << i << " " << splice_rsd[ i ] << std::endl;
	// }

	// create subpose of chain TM regions from the PDB into a new pose
	core::pose::Pose subpose = core::pose::Pose();
	pdbslice( subpose, pose, splice_rsd );

	// subpose.dump_pdb("subpose.pdb");
	// pose.dump_pdb("pose.pdb");

	// compute COM of the newly created pose
	core::Vector com( get_center_of_mass( subpose ) );

	return com;

} // chain TM COM

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////

/// @brief Residue closest to pose transmembrane center-of-mass
/// @details Gets the coordinates of the residue closest to TM span center-of-mass
///   This only looks at the span start and end residues for
///   calculation of the TM span COM, this should be faster than the
///   real thing though
core::Size rsd_closest_to_pose_tm_com( core::pose::Pose const & pose ) {

	// get pose TM com
	core::Vector com( pose_tm_com( pose ) );

	// get closest residue
	int nearest( return_nearest_residue( pose, 1, nres_protein( pose ), com ) );

	return static_cast< core::Size >( nearest );

} // residue closest to pose TM COM

/////////////////////////////////////////

/// @brief Residue closest to chain center-of-mass
/// @details Gets the residue number closest to the chain center-of-mass
core::Size rsd_closest_to_chain_com( core::pose::Pose const & pose, int chainid ) {

	// check that the chain isn't the membrane residue
	if ( chainid == pose.chain( pose.conformation().membrane_info()->membrane_rsd_num() ) ) {
		utility_exit_with_message("You are trying to compute the center-of-mass for the membrane residue as a chain. Choose a different one...");
	}

	// find start and end residue number for chain in question
	int start( 1 );
	int end( nres_protein( pose ) );

	for ( core::Size i = 1; i <= nres_protein( pose ); ++i ) {

		// find starting position
		if ( i > 1 && pose.chain( i ) == chainid && pose.chain( i-1 ) != chainid ) {
			start = i;
		}

		// find end position
		if ( i < nres_protein( pose ) && pose.chain( i ) == chainid && pose.chain( i+1 ) != chainid ) {
			end = i;
		}
	}

	TR << "start: " << start << " end: " << end << std::endl;

	// compute chain COM
	core::Size rsd_id = static_cast< core::Size >( residue_center_of_mass( pose, start, end ) );

	return rsd_id;

} // residue closest to chain COM

////////////////////////////////////////////////////////////////////////////////

/// @brief Residue closest to chain TM center-of-mass
/// @details Gets the residue number closest to the chain TM center-of-mass
core::Size rsd_closest_to_chain_tm_com( core::pose::Pose const & pose, int chainid ) {

	// check for ligand:
	// find out the residue number belonging to ligand, the entire chain is ligand
	core::Size chain_end_rsd( chain_end_res( pose, static_cast< core::Size >( chainid ) ) );

	// if residue is ligand, then return ligand residue number
	if ( ! pose.residue( chain_end_rsd ).is_protein() ) {
		TR << "chain end rsd " << chain_end_rsd << " is not protein" << std::endl;
		return chain_end_rsd;
	}

	// get membrane residue number and chain for checking
	core::Size memrsd = pose.conformation().membrane_info()->membrane_rsd_num();
	int mem_chain = pose.residue( memrsd ).chain();

	if ( chainid == mem_chain ) {
		utility_exit_with_message( "You are trying to get the rsd closest to chain TM COM, but you are looking at the membrane residue. Quitting." );
	}

	// compute chain TM com
	TR << "chainid " << chainid << " chain " << get_chain_from_chain_id( static_cast<core::Size >( chainid ), pose ) << std::endl;
	core::Vector com = chain_tm_com( pose, chainid );

	TR << "com " << com.to_string() << std::endl;

	// go through each residue in the pose and find CA atom that is closest
	core::Real min_dist = 999999.0;
	core::Size min_rsd = 0;
	for ( core::Size i = 1; i <= nres_protein( pose ); ++i ) {

		// get CA position of this residue and compute distance to COM
		core::Vector ca_pos = pose.residue( i ).atom( "CA" ).xyz() ;
		core::Real distance = ( ca_pos - com ).length();

		// get chainid and make sure we only pick a residue within that chain
		int rsd_chainid = pose.residue( i ).chain();

		if ( distance < min_dist && rsd_chainid == chainid ) {
			min_dist = distance;
			min_rsd = i;
			//   TR << "min dist " << min_dist << " min rsd " << i << std::endl;
		}
	}

	return min_rsd;

} // rsd_closest_to_chain_tm_com

////////////////////////////////////////////////////////////////////////////////

/// @brief Check reasonable range of vector
void check_vector( core::Vector const vector ) {

	TR << "Checking vector " << std::endl;

	// warn if vector is origin
	if ( vector.to_string() == "(0, 0, 0)" ) {
		TR << "WARNING: your vector is (0, 0, 0)!" << std::endl;
	}

	// Fail if vector is out of range
	if ( vector.x() < -1000 || vector.x() > 1000 ||
			vector.y() < -1000 || vector.y() > 1000 ||
			vector.z() < -1000 || vector.z() > 1000 ) {

		throw new core::conformation::membrane::EXCN_Illegal_Arguments("Unreasonable range for center or normal! Check your input vectors!");
	}
}// check_vector

/// @brief Normalize normal vector to length 15 for visualization
void membrane_normal_to_length_15( core::pose::Pose & pose ){

	using namespace core;

	// get center and normal
	Vector center = pose.conformation().membrane_info()->membrane_center();
	Vector normal = pose.conformation().membrane_info()->membrane_normal();

	// normalize normal vector
	normal.normalize( 15 );

	// Update membrane position with new coords
	pose.conformation().update_membrane_position( center, normal );
}

/// @brief Calculates translation axis lying in the membrane (= projection axis
///   between embedding centers)
core::Vector const membrane_axis( core::pose::Pose & pose, int jumpnum )
{
	using namespace core::pose;
	using namespace numeric;
	using numeric::cross;
	using namespace protocols::membrane::geometry;

	// get embedding between partners
	EmbeddingDef emb_up, emb_down;
	update_partner_embeddings( pose, static_cast< core::Size > ( jumpnum ), emb_up, emb_down );

	// temporary axis is axis between partner embedding centers
	// this doesn't need to be directly in the membrane plane.
	// since we need this to be exactly in the membrane plane (since moving
	// apart 100A or more can easily move things out of the membrane), we will
	// compute the projection axis of the tmp_axis onto the membrane plane
	core::Vector tmp_axis = emb_up.center() - emb_down.center();

	// get membrane normal
	core::Vector mem_normal = pose.conformation().membrane_info()->membrane_normal();

	// compute axis orthogonal to both tmp_axis and membrane normal
	core::Vector ortho = cross( mem_normal, tmp_axis );

	// compute axis orthogonal to ortho and membrane normal
	// core::Vector trans_axis = numeric::cross( ortho, mem_normal );
	return cross( ortho, mem_normal );

}// membrane axis

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
	core::pose::Pose const & pose,     // full pose
	core::Size const jumpnum,     // jump number to split on
	core::conformation::membrane::SpanningTopology const & topo,  // topology to split
	core::pose::Pose & pose_up,      // upstream partner after pose splitting
	core::pose::Pose & pose_down,     // downstream partner after pose splitting
	core::conformation::membrane::SpanningTopology & topo_up,   // topology of upstream pose
	core::conformation::membrane::SpanningTopology & topo_down  // topology of downstream pose
) {
	// can't split pose by membrane jump, partition_pose_by_jump function will fail
	if ( jumpnum == static_cast< core::Size > ( pose.conformation().membrane_info()->membrane_jump() ) ) {
		utility_exit_with_message("Cannot split pose by membrane jump! Quitting...");
	}

	// split pose at jump, both new poses have residue numbering starting at 1!!!
	partition_pose_by_jump( pose, jumpnum, pose_up, pose_down );

	// go through TMspans
	for ( core::Size i = 1; i <= topo.nspans(); ++i ) {

		// get start and end residues
		core::Size start = topo.span( i )->start();
		core::Size end = topo.span( i )->end();

		// if span is not in upstream partner, add to downstream topology
		if ( start > pose_up.total_residue() ) {

			core::Size new_start = start - nres_protein( pose_up );
			core::Size new_end = end - nres_protein( pose_up );

			// add to downstream topology
			topo_down.add_span( new_start, new_end, 0 );
		} else {
			// else add to upstream topology
			topo_up.add_span( start, end, 0 );
		}
	}

	TR << "topology of upward partner" << std::endl;
	topo_up.show();
	TR << "topology of downward partner" << std::endl;
	topo_down.show();
}// split topology by jump


/// @brief Splits the SpanningTopology object into two objects, depending on
/// given jump number
/// @details This doesn't shift the topology to start at 1 for each partner, it
/// remains exactly the same as it would be for the complete pose, just split
/// BEWARE: this does not work for splitting topology by spans! It only works
/// chainwise
void split_topology_by_jump_noshift(
	core::pose::Pose const & pose,  // full pose
	core::Size const jumpnum,    // jump number to split on
	core::conformation::membrane::SpanningTopologyOP topo, // topology to split
	core::conformation::membrane::SpanningTopologyOP topo_up,  // topology of upstream pose
	core::conformation::membrane::SpanningTopologyOP topo_down // topology of downstream pose
) {
	// can't split pose by membrane jump, partition_pose_by_jump function will fail
	if ( jumpnum == static_cast< core::Size > ( pose.conformation().membrane_info()->membrane_jump() ) ) {
		utility_exit_with_message("Cannot split pose by membrane jump! Quitting...");
	}

	// MAKING ASSUMPTION THAT DOWNSTREAM PARTNER IS A SINGLE CHAIN!!!
	// alternative implementation without making that assumption:
	// = implement a function in the FoldTree that gives you all residues downstream
	//   of a jump
	// = go through entire pose and write a vector1<bool> per residue saying true
	//   for downstream partner and false for upstream

	// find chain of downstream residue number
	core::Size res_downstream = static_cast< core::Size > ( pose.fold_tree().downstream_jump_residue( jumpnum ) );
	int chain_downstream = pose.chain( res_downstream );

	// go through TMspans
	for ( core::Size i = 1; i <= topo->nspans(); ++i ) {

		// get start and end residues
		core::Size start = topo->span( i )->start();
		core::Size end = topo->span( i )->end();

		// if span is not in upstream partner, add to downstream topology
		if ( pose.chain( start ) == chain_downstream ) {
			topo_down->add_span( start, end, 0 );
		} else {
			// else add to upstream topology
			topo_up->add_span( start, end, 0 );
		}
	}

}// split topology by jump, no shift in topology objects

/// @brief Split topology by chain
/// @details Split topology by chain and give vector of topology objects
utility::vector1< core::conformation::membrane::SpanningTopologyOP >
split_topology_by_chain_noshift(
	core::pose::Pose const & pose,
	core::conformation::membrane::SpanningTopologyOP const topo
) {

	using namespace core::conformation::membrane;

	// initialize variables
	utility::vector1< SpanningTopologyOP > topos;

	// go through each chain
	for ( int c = 1; c <= pose.chain( nres_protein( pose ) ); ++c ) {

		// create new topology object
		SpanningTopologyOP chain_topo( new SpanningTopology() );

		// go through total spanning topology
		for ( core::Size s = 1; s <= topo->nspans(); ++s ) {

			// get start and end
			core::Size start = topo->span( s )->start();
			core::Size end = topo->span( s )->end();

			// get chain
			int chain = pose.chain( start );

			// if span is part of the chain, add it to chain topology
			if ( chain == c ) {
				chain_topo->add_span( start, end );
			}

		} // spans

		// add chain topologies to total topology
		topos.push_back( chain_topo );

	} // chains

	return topos;

} // split topology by chain


/////////////////////////////////////////////
// Methods for reading center/normal in IO //
/////////////////////////////////////////////

/// @brief Read in a user provided center/normal pair from RosettaScripts
/// @details Given an XML tag from a RosettaScript read in a center & normal
/// option into two xyzVector objects. This method is intended to reduce duplication
/// accross membrane framework movers that use the same tricks for vectors.
/// Takes two Vector references and a Tag&
void
read_center_normal_from_tag( core::Vector & center, core::Vector & normal, utility::tag::TagCOP tag ) {

	using namespace core;

	// Read in membrane center & normal
	if ( tag->hasOption( "center" ) ) {
		std::string cen = tag->getOption< std::string >( "center" );
		utility::vector1< std::string > str_cen = utility::string_split_multi_delim( cen, ":,'`~+*&|;." );

		if ( str_cen.size() != 3 ) {
			utility_exit_with_message( "Cannot read in xyz center vector from string - incorrect length!" );
		} else {
			center.x() = std::atof( str_cen[1].c_str() );
			center.y() = std::atof( str_cen[2].c_str() );
			center.z() = std::atof( str_cen[3].c_str() );
		}
	}

	if ( tag->hasOption( "normal" ) ) {
		std::string norm = tag->getOption< std::string >( "normal" );
		utility::vector1< std::string > str_norm = utility::string_split_multi_delim( norm, ":,'`~+*&|;." );

		if ( str_norm.size() != 3 ) {
			utility_exit_with_message( "Cannot read in xyz center vector from string - incorrect length!" );
		} else {
			normal.x() = std::atof( str_norm[1].c_str() );
			normal.y() = std::atof( str_norm[2].c_str() );
			normal.z() = std::atof( str_norm[3].c_str() );
		}
	}
}

/// @brief Read in a user provided center/normal pair from the commandline, safetly
/// @details Read the membrane setup center & normal options form the command line
/// from mp:setup:center and mp:setup:normal. Intended to reduce IO code duplication
/// in membrane framework movers.
void
read_center_normal_from_cmd( core::Vector & center, core::Vector & normal ) {

	using namespace basic::options;

	// Read in Center Parameter
	if ( option[ OptionKeys::mp::setup::center ].user() ) {
		if ( option[ OptionKeys::mp::setup::center ]().size() == 3 ) {
			center.x() = option[ OptionKeys::mp::setup::center ]()[1];
			center.y() = option[ OptionKeys::mp::setup::center ]()[2];
			center.z() = option[ OptionKeys::mp::setup::center ]()[3];
		} else {
			utility_exit_with_message( "Center xyz vector must have three components! Option has either too many or too few arguments!" );
		}
	}

	// Read in Normal Parameter
	if ( option[ OptionKeys::mp::setup::normal ].user() ) {
		if ( option[ OptionKeys::mp::setup::normal ]().size() == 3 ) {
			normal.x() = option[ OptionKeys::mp::setup::normal ]()[1];
			normal.y() = option[ OptionKeys::mp::setup::normal ]()[2];
			normal.z() = option[ OptionKeys::mp::setup::normal ]()[3];
		} else {
			utility_exit_with_message( "Normal xyz vector must have three components! Option has either too many or too few arguments" );
		}
	}

}

} // membrane
} // protocols
