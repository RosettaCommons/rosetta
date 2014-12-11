// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		protocols/membrane/geometry/util.cc
///
/// @brief 		Utility methods for membrane framework
/// @details 	Utility methods include determining center of mass (moved down in the tree)
///				and adjusting normal parameters for visualization.
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_geometry_util_cc
#define INCLUDED_protocols_membrane_geometry_util_cc

// Unit Headers
#include <protocols/membrane/geometry/util.hh>

// Project Headers
#include <core/conformation/membrane/Exceptions.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/Embedding.hh>

// Package Headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility Headers
#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

#include <basic/Tracer.hh>

// C++ Headers
#include <algorithm>
#include <string>
#include <cstdlib>
#include <cmath>

static thread_local basic::Tracer TR( "protocols.membrane.geometry.util" );

using basic::Error;
using basic::Warning;

using namespace core;
using namespace core::conformation::membrane;

namespace protocols {
namespace membrane {
namespace geometry {


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
		   pose::Pose const & pose,
		   core::SSize const start,
		   core::SSize const stop
		   )
{
	Vector center( 0.0 );
	for ( core::SSize i = start; i <= stop; ++i ) {
		if( !pose.residue( i ).is_protein() ) {
			Vector ca_pos( pose.residue( i ).nbr_atom_xyz() );
			center += ca_pos;
		} else {
			Vector ca_pos( pose.residue( i ).atom( "CA" ).xyz() );
			center += ca_pos;
		}
	}
	center /= ( stop - start + 1 );

	return center;
}


/// @brief      Residue Center of Mass
/// @details    Calcualte the center of mass of a pose.
///
/// @author     Monica Berrondo, Modified by Javier Castellanos and Rebecca Alford
core::SSize
residue_center_of_mass(
				   pose::Pose const & pose,
				   core::SSize const start,
				   core::SSize const stop
				   )
{
	Vector center = center_of_mass( pose, start, stop );
	return return_nearest_residue( pose, start, stop, center );
}

/// @brief      Return nearest residue
/// @details    Find the residue nearest some position passed in (normally a center of mass)
///
/// @author     Monica Berrondo, Modified by Javier Castellanos and Rebecca Alford
core::SSize
return_nearest_residue(
				   pose::Pose const & pose,
				   core::SSize const begin,
				   core::SSize const end,
				   Vector center
				   )
{
	Real min_dist = 9999.9;
	core::SSize res = 0;
	for ( core::SSize i=begin; i<=end; ++i )
	{
		Vector ca_pos;
		if ( !pose.residue( i ).is_protein() ){
			ca_pos = pose.residue( i ).nbr_atom_xyz();
		} else {
			ca_pos = pose.residue( i ).atom( "CA" ).xyz() ;
		}
		
		ca_pos -= center;
		Real tmp_dist( ca_pos.length_squared() );

		if ( tmp_dist < min_dist ) {
			res = i;
			min_dist = tmp_dist;
		}
	}
	return res;
}

/// @brief		Get z-coord and chainID
/// @details	Helper function that creates input for SpanningTopology
///				which is not built at the time the Pose is built
///				returns a pair of vectors:
///				vector1 is z-coord of CA atoms of the pose
///				vector2 is chainID of CA atoms of the pose
std::pair< utility::vector1< Real >, utility::vector1< Real > >
get_chain_and_z( pose::PoseOP pose ) {
	
	TR.Debug << "get_pose_info" << std::endl;
	using namespace core::pose;
	
	// initialize variables
	utility::vector1< Real > z_coord;
	utility::vector1< Size > chain_info;
	
	// loop over residues, get chain info and z_coord
	for ( Size i = 1; i <= pose->total_residue(); ++i ){
		
		// get info
		z_coord.push_back( static_cast< Real >(pose->residue(i).atom(2).xyz().z()) );
		chain_info.push_back( pose->chain(i) );
	}
	
	// put the data in a pair
	std::pair< utility::vector1< Real >, utility::vector1< Size > > pose_info( z_coord, chain_info );
	
	return pose_info;
	
} // get chain and z from pose

////////////////////////////////////////////////////////////////////////////////
/// @brief Compute Membrane Center/Normal from Membrane Spanning
/// topology
void compute_structure_based_membrane_position(
										pose::Pose & pose,
										Vector & center,
										Vector & normal
										) {
	// create EmbeddingDef to return
	EmbeddingDefOP embed = compute_structure_based_membrane_position( pose );

	// set new center and normal
	center = embed->center();
	normal = embed->normal();
	
}// compute structure-based position

////////////////////////////////////////////////////////////////////////////////
/// @brief Compute Membrane Center/Normal from Membrane Spanning
/// topology
EmbeddingDefOP compute_structure_based_membrane_position( pose::Pose & pose ){

	using namespace protocols::membrane::geometry;
	using namespace core::conformation::membrane;
	
	// get SpanningTopology
	SpanningTopologyOP topo = pose.conformation().membrane_info()->spanning_topology();
	
	// get create Embedding object
	EmbeddingOP embeddings( new Embedding( topo, pose ) );
	
	// return total embedding
	return embeddings->total_embed();
			
}// compute structure based membrane position

/// @brief Check reasonable range of vector
void check_vector( core::Vector const vector ) {

	TR << "Checking vector " << std::endl;

	// warn if vector is origin
	if ( vector.to_string() == "(0, 0, 0)"){
	TR << "WARNING: your vector is (0, 0, 0)!" << std::endl;
	}

	// Fail if vector is out of range
	if ( vector.x() < -1000 || vector.x() > 1000 ||
	vector.y() < -1000 || vector.y() > 1000 ||
	vector.z() < -1000 || vector.z() > 1000 ) {

		throw new conformation::membrane::EXCN_Illegal_Arguments("Unreasonable range for center or normal! Check your input vectors!");
	}
}// check_vector

/// @brief Average EmbeddingDefs as they are without vector inversion accounting for topology
/// @details Get average center and normal from a vector of EmbeddingDefs
EmbeddingDefOP average_embeddings( utility::vector1< EmbeddingDefOP > parts ) {

    // Initialize vars
    core::Vector center(0, 0, 0);
    core::Vector normal(0, 0, 0);
    
	// Compute resulting center and normal
    for ( Size i = 1 ; i <= parts.size(); ++i ) {
        center += parts[i]->center();
		normal += parts[i]->normal();
    }
	
    center /= parts.size();
	normal.normalize( 15 );
    
	// Create new embedding setup and return it
    EmbeddingDefOP embedding( new EmbeddingDef( center, normal ) );
    return embedding;

}// average embeddings

/// @brief Average EmbeddingDefs after first inverting some vectors accounting for topology
/// @details Get average center and normal from a vector of EmbeddingDefs
EmbeddingDefOP average_antiparallel_embeddings( utility::vector1< EmbeddingDefOP > parts ) {

	// Initialize vars
	core::Vector center(0, 0, 0);
	core::Vector normal(0, 0, 0);

	// embedding of first span
	Vector const center1 = parts[1]->center();
	Vector const normal1 = parts[1]->normal();

	// Compute resulting center and normal
	for ( Size i = 1 ; i <= parts.size(); ++i ) {

		TR << "center: " << parts[i]->center().to_string() << "normal: " << parts[i]->normal().to_string() << std::endl;

		// calculate new center
		center += parts[i]->center();
		
		// calculate points for angle calculation
		Vector p1 = center1 + normal1;
		Vector p  = center1 + parts[i]->normal();
		
		// calculate  angle between normals of first object and this one
		Real angle( numeric::angle_degrees( p1, center1, p ) );
		TR << "angle: " << angle << std::endl;
		
		// check if angle of normal is < 100 degrees to first normal
		// if yes, then add to normal, if no add inverted vector
		if ( angle > -100 && angle < 100 ) {
			normal += parts[i]->normal();
		}
		else {
			normal -= parts[i]->normal();
		}
	}

	center /= parts.size();
	normal.normalize( 15 );
	
	// Create new embedding setup and return it
	EmbeddingDefOP embedding( new EmbeddingDef( center, normal ) );
	return embedding;
	
}// average antiparallel embeddings

/// @brief Normalize normal vector to length 15 for visualization
void membrane_normal_to_length_15( pose::Pose & pose ){

	// get center and normal
	Vector center = pose.conformation().membrane_info()->membrane_center();
	Vector normal = pose.conformation().membrane_info()->membrane_normal();

	// normalize normal vector
	normal.normalize( 15 );

	// Update membrane position with new coords
	pose.conformation().update_membrane_position( center, normal );
}

/// @brief Set membrane residue to root of foldtree
/// @details Requires MembraneInfo to be constructed beforehand;
///			 use AddMembraneMover to do that
void reorder_membrane_foldtree( pose::Pose & pose ) {

	// get foldtree from pose
	core::kinematics::FoldTree foldtree = pose.fold_tree();

	// reorder foldtree
	foldtree.reorder( pose.conformation().membrane_info()->membrane_rsd_num() );

	// set foldtree in pose
	pose.fold_tree( foldtree );
}

} // geometry
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_geometry_util_cc

