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
#include <protocols/membrane/geometry/EmbeddingDef.hh>

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

static thread_local basic::Tracer TR( "protocols.membrane.util" ); 

using basic::Error;
using basic::Warning;

namespace protocols {
namespace membrane {
    
using namespace numeric;
using namespace core;
using namespace core::scoring;
using namespace core::pose;
using namespace core::conformation::membrane;
using namespace protocols::membrane::geometry;

/////////////////////////////////////////////////////////////////////////
// Methods for calculating rmsds between protein transmembrane regions //
/////////////////////////////////////////////////////////////////////////

/// @brief Compute backbone RMSD between TM regions - don't superimpose
/// @details Calculate the rmsd between backbone atoms (N, CB, CA, O)
/// in the transmembrane regions, as defined by the spanning topology
/// object Do not superimpose the poses. Takes a native pose and
/// current pose
core::Real
mem_bb_rmsd_no_super( Pose & native_pose, Pose & pose ) {
    
    // Check that pose is actually a membrane pose
    if (! pose.conformation().is_membrane()) {
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
mem_all_atom_rmsd_no_super( Pose & native_pose, Pose & pose ) {
    
    // Check that pose is actually a membrane pose
    if (! pose.conformation().is_membrane()) {
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
    
    return rmsd_no_super_subset( native_pose, pose, tm_regions, is_heavyatom );
    
}

/// @brief Compute backbone RMSD between TM regions - do superimpose
/// @details Calculate the rmsd between backbone atoms (N, CB, CA, O)
/// in the transmembrane regions, as defined by the spanning topology
/// object Superimpose the poses. Takes a native pose and current pose
core::Real
mem_bb_rmsd_with_super( Pose & native_pose, Pose & pose ) {
    
    // Check that pose is actually a membrane pose
    if (! pose.conformation().is_membrane()) {
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
mem_all_atom_rmsd_with_super( Pose & native_pose, Pose & pose )  {
    
    // Check that pose is actually a membrane pose
    if (! pose.conformation().is_membrane() ) {
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
calc_helix_tilt_angle( Pose & pose, core::Size span_no ) {
    
    // Get membrane normal from membrane info
    Vector normal( pose.conformation().membrane_info()->membrane_normal() );
	
	// Calculate an axis representing the "axis of the helix"
    Vector helix_axis( calc_helix_axis( pose, span_no ) );
    
    // Normalize both vectors (just incase
    normal.normalize();
    helix_axis.normalize();
    
    // Check that both vectors have a positive orientation
    if ( normal.length() > 0 ) normal.negate();
    if ( helix_axis.length() > 0 ) helix_axis.negate();
    
    // Calculate the angle between the normal & helix
    core::Real tilt_angle( angle_of( normal, helix_axis ) );
    return numeric::conversions::degrees( tilt_angle );
}
    
/// @brief Determine the axis used to define a single TM Helix
/// @details Using the COM of the helix start & end position, calculate a helix
/// describing its geometry relative to the memrbane normal. Takes a pose &
/// span number. Not a good approx for helices with kinks.
core::Vector
calc_helix_axis( Pose & pose, core::Size span_no ) {
    
    // Get the span from the pose
    SpanOP helix_span( pose.conformation().membrane_info()->spanning_topology()->span( span_no ) );
	
	// Check the size of the span is sufficient for the this calculation
	if ( helix_span->end() - helix_span->start() < 6 ) {
		utility_exit_with_message( "Transmembrane span is too small to calculate a helix axis using the center of masses - less than six residues" );
	}
	
	// Check the size of the span relative to the size of the protein is appropriate for
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
    
////////////////////////////////////////////////////////////////
// Safety checks & convenience methods for membrane foldtrees //
////////////////////////////////////////////////////////////////

/// @brief Determine whether the membrane is modeled as fixed
/// @details Based on the setup of the foldtree, determined whether
/// the membrane is currently fixed, meaning it is setup as the
/// root in the FoldTree and has no upstream children. Takes a pose.
bool
is_membrane_fixed( Pose & pose ) {

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
is_membrane_moveable_by_itself( Pose & pose ) {
    
    using namespace core::kinematics;
	
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
        for ( kinematics::FoldTree::const_iterator it = current_ft.begin(), it_end = current_ft.end(); it != it_end; ++it ) {
            
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
void reorder_membrane_foldtree( pose::Pose & pose ) {
    
    // get foldtree from pose
    core::kinematics::FoldTree foldtree = pose.fold_tree();
    
    // reorder foldtree
    foldtree.reorder( pose.conformation().membrane_info()->membrane_rsd_num() );
    
    // set foldtree in pose
    pose.fold_tree( foldtree );
}

///////////////////////////////////////////////////////////
// Utilities for accessing dssp, z coords and chain info //
///////////////////////////////////////////////////////////

/// @brief Grab the z-coordinates and chainIDs from the entire pose
/// @details From the pose, grab all of the z_coords of CA atoms and
/// chain IDs, currently used for spanning topology construction.
/// Returns a std::pair of two vectors: the first a vector1 of z
/// coordinates and the second a vector1 of chainIDs for CA atoms
std::pair< utility::vector1< core::Real >, utility::vector1< core::Real > >
get_chain_and_z( pose::Pose const & pose ) {
    
    TR.Debug << "get_pose_info" << std::endl;
    using namespace core::pose;
    
    // initialize variables
    utility::vector1< core::Real > z_coord;
    utility::vector1< core::Size > chain_info;
    
    // loop over residues, get chain info and z_coord
    for ( core::Size i = 1; i <= pose.total_residue(); ++i ){
        
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
utility::vector1< char > get_secstruct( pose::Pose & pose ) {
    
    TR.Debug << "get_secstruct" << std::endl;
    using namespace core::pose;
    
    // set secondary structure in pose with DSSP
    protocols::moves::DsspMover dssp = protocols::moves::DsspMover();
    dssp.apply( pose );
    
    // initialize variable
    utility::vector1< char > secstruct;
    
    // loop over residues, get chain info and z_coord
    for ( core::Size i = 1; i <= pose.total_residue(); ++i ){
        
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
       pose::Pose const & pose,
       SpanningTopology const & topology,
       Vector & center,
       Vector & normal
       ) {
    
    // create EmbeddingDef to return
    EmbeddingDefOP embed = compute_structure_based_embedding( pose, topology );
    
    // set new center and normal
    center = embed->center();
    normal = embed->normal();
    
}// compute structure-based position

    
/// @brief Compute Membrane Center/Normal from Membrane Spanning
/// topology, uses topology from MembraneInfo
void compute_structure_based_embedding(
    pose::Pose const & pose,
    Vector & center,
    Vector & normal
    ) {
    
    // get topology from MembraneInfo
    SpanningTopology topo( *pose.conformation().membrane_info()->spanning_topology() );
    
    // create EmbeddingDef to return
    compute_structure_based_embedding( pose, topo, center, normal );
    
}// compute structure-based position

/// @brief Compute Membrane Center/Normal from Membrane Spanning
/// topology
EmbeddingDefOP compute_structure_based_embedding( pose::Pose const & pose, SpanningTopology const & topo ){
    
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
EmbeddingDefOP compute_structure_based_embedding( pose::Pose const & pose ){
    
    // get topology from MembraneInfo
    SpanningTopology topo( *pose.conformation().membrane_info()->spanning_topology() );
    
    return compute_structure_based_embedding( pose, topo );
    
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
EmbeddingDefOP average_embeddings( utility::vector1< EmbeddingDefOP > const parts ) {
    
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
EmbeddingDefOP average_antiparallel_embeddings( utility::vector1< EmbeddingDefOP > const parts ) {
    
    // Initialize vars
    core::Vector center(0, 0, 0);
    core::Vector normal(0, 0, 0);
    
    // embedding of first span
    Vector const center1 = parts[1]->center();
    Vector const normal1 = parts[1]->normal();
    
    // Compute resulting center and normal
    for ( core::Size i = 1 ; i <= parts.size(); ++i ) {
        
        TR << "center: " << parts[i]->center().to_string() << "normal: " << parts[i]->normal().to_string() << std::endl;
        
        // calculate new center
        center += parts[i]->center();
        
        // calculate points for angle calculation
        Vector p1 = center1 + normal1;
        Vector p  = center1 + parts[i]->normal();
        
        // calculate  angle between normals of first object and this one
        core::Real angle( numeric::angle_degrees( p1, center1, p ) );
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
    normal.normalize();
    
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

/// @brief Calculates translation axis lying in the membrane (= projection axis
///			between embedding centers)
core::Vector const membrane_axis( pose::Pose & pose, int jumpnum )
{
    using namespace core::pose;
    using namespace numeric;
    using numeric::cross;
    
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
///	given jump number
/// @details This is useful for calculating an embedding for parts of the
///	structure: this can now easily be accomplished by creating two empty topology
/// objects, call this function, and then use both topology objects and subposes
/// to call compute_structure_based_membrane_embedding
///	BEWARE: this does not work for splitting topology by spans! It only works
///	chainwise
void split_topology_by_jump(
    Pose const & pose,					// full pose
    core::Size const jumpnum,					// jump number to split on
    SpanningTopology const & topo,		// topology to split
    Pose & pose_up,						// upstream partner after pose splitting
    Pose & pose_down,					// downstream partner after pose splitting
    SpanningTopology & topo_up,			// topology of upstream pose
    SpanningTopology & topo_down		// topology of downstream pose
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
        }
        
        // else add to upstream topology
        else {
            topo_up.add_span( start, end, 0 );
        }
    }
    
    TR << "topology of upward partner" << std::endl;
    topo_up.show();
    TR << "topology of downward partner" << std::endl;
    topo_down.show();
}// split topology by jump

    
/// @brief Splits the SpanningTopology object into two objects, depending on
///	given jump number
/// @details This doesn't shift the topology to start at 1 for each partner, it
///	remains exactly the same as it would be for the complete pose, just split
///	BEWARE: this does not work for splitting topology by spans! It only works
///	chainwise
void split_topology_by_jump_noshift(
    pose::Pose const & pose,		// full pose
    core::Size const jumpnum,				// jump number to split on
    SpanningTopologyOP topo,	// topology to split
    SpanningTopologyOP topo_up,		// topology of upstream pose
    SpanningTopologyOP topo_down	// topology of downstream pose
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
        }
        
        // else add to upstream topology
        else {
            topo_up->add_span( start, end, 0 );
        }
    }
    
}// split topology by jump, no shift in topology objects

/// @brief Update embedding of the partners after a move
/// @details Requires the jump number between the partners, the topology will
///				be taken from MembraneInfo and will be split accordingly; up and
///				down means upstream and downstream
void update_partner_embeddings( pose::Pose const & pose, core::Size const jumpnum, EmbeddingDef & emb_up, EmbeddingDef & emb_down ) {
    
    // SpanningTopology objects
    SpanningTopologyOP topo = pose.conformation().membrane_info()->spanning_topology();
    SpanningTopologyOP topo_up( new SpanningTopology() );	// upstream partner
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
    
/////////////////////////////////////////////
// Methods for reading center/normal in IO //
/////////////////////////////////////////////

/// @brief Read in a user provided center/normal pair from RosettaScripts
/// @details Given an XML tag from a RosettaScript read in a center & normal
/// option into two xyzVector objects. This method is intended to reduce duplication
/// accross membrane framework movers that use the same tricks for vectors.
/// Takes two Vector references and a Tag&
void
read_center_normal_from_tag( Vector & center, Vector & normal, utility::tag::TagCOP tag ) {
    
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
read_center_normal_from_cmd( Vector & center, Vector & normal ) {
    
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
