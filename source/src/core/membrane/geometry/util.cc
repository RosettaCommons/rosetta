// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/geometry/util.cc
///
/// @brief 		Utility methods for defining membranes and membrane embeddings
/// @details 	Helps to check for internal errors, bounds, and object equality
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_geometry_util_cc
#define INCLUDED_core_membrane_geometry_util_cc

// Unit Headers
#include <core/membrane/geometry/util.hh>

// Project Headers
#include <core/membrane/properties/SpanningTopology.hh>
#include <core/membrane/util/Exceptions.hh>

// Package Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>

#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

#include <basic/Tracer.hh>

// C++ Headers
#include <algorithm>
#include <string>
#include <cstdlib>
#include <cmath>

static basic::Tracer TR("core.membrane.geometry.util");

using basic::Error;
using basic::Warning;

using namespace core::membrane::properties;

namespace core {
namespace membrane {
namespace geometry {

    /// @brief      Virtual residue Equals
    /// @details    Custom equality method - Checks two virtual atoms are equal in typesets
    ///             and cartesian coordinates
    ///
    /// Precondition: (Checked) Both residues are virtual residues
    ///
    /// @param 	rsd1
    ///				first atom to investigate
    /// @param 	rsd2
    ///				second residue to investigate
    /// @param fullatom
    ///				specifies if the atom typesef of the pose is fullatom
    ///
    ///	@return	bool
    bool virtual_rsd_equal( core::conformation::ResidueOP rsd1, core::conformation::ResidueOP rsd2 ) {

        // The first two conditions check that each residue is a virtual residue
        // which allows us to surpass additiona typeset checking. These are fatal cases
        // indicating misuse of the method and will cause Rosetta to fail.

        // Check that the first residue is virtual
        if ( rsd1->aa() != core::chemical::aa_vrt ) { return false; }
        if ( rsd2->aa() != core::chemical::aa_vrt ) { return false; }

        // Checks virtual atom coordinates are equal
        if ( rsd1->atom(2).xyz().x() != rsd2->atom(2).xyz().x() ) { return false; }
        if ( rsd1->atom(2).xyz().y() != rsd2->atom(2).xyz().y() ) { return false; }
        if ( rsd1->atom(2).xyz().z() != rsd2->atom(2).xyz().z() ) { return false; }
        
        // All of the conditions have been met!
        return true;
    }

    /// @brief      Get Residue Depth in Membrane
    /// @details    Calculate the depth of a residue with respect to membrane players
    ///
    /// @param      normal
    ///                 provided normal vector to membrane
    /// @param      center
    ///                 provided center point for embedding for membrane
    core::Real
    get_mpDepth( core::Vector normal, core::Vector center, core::conformation::Residue rsd ) {
        
        // Get membrane xyz
        Vector const & xyz( rsd.atom( 2 ).xyz() );
        core::Real depth = dot( xyz-center, normal )+30;
        return depth;
    }
    
    /// @brief      Check Membrane Spanning
    /// @details    Check that caucluated membrane spanning respects new
    ///             normal and center definitions
    ///
    /// @throws     <none>
    /// @note       Needs refactoring!!!
	bool
    check_spanning(
                   core::pose::Pose const & pose,
                   core::Vector const & normal,
                   core::Vector const & center,
                   SpanningTopologyOP topology
                   ) {
        
        // Loop Through the Pose
        for( Size i = 1; i <= topology->total_tmhelix()-1 ; ++i ) {
            
            // If scoring allowed at position, continue
            if ( ! topology->allow_tmh_scoring()[i] ) continue;
            
            Vector const & start_i( pose.residue( topology->span()(i, 1) ).atom( 2 ).xyz());
            bool start_i_side=(dot(start_i-center,normal) > 0);
            bool span_check=false;
            
            // Loop through iteratively
            for( Size j = i+1; j <= topology->total_tmhelix(); ++j) {
                
                // Check scoring is allowed at given position
                if( !topology->allow_tmh_scoring()[j] || span_check ) continue;
                span_check=true;
                
                Vector const & start_j( pose.residue( topology->span()(j, 2) ).atom( 2 ).xyz());
                bool start_j_side=(dot(start_j-center,normal) > 0);
                bool coord_para=(start_i_side==start_j_side);
                
                if ( topology->helix_id()[i]-topology->helix_id()[i] % 2 == 0 ) {
                    
                    if(!(coord_para)) { return false; }
                    
                } else {
                    
                    if(coord_para) { return false; }
                }
            }
        }
        return true;
    }
    
    ///// Embedding Calculation Utility Mehtods //////////
    
    /// @brief    Sum 2 Numeric XYZ Reals
    /// @details  Calculate the vector sum of 2 nymeric xyz vectors
    ///           utility function for the two methods below - maintains precision for core::Real
    numeric::xyzVector< core::Real >
    xyz_sum( numeric::xyzVector< core::Real > & a, numeric::xyzVector< core::Real > & b) {
        
        numeric::xyzVector< core::Real > sum;
        core::Real sum_x = a.x() + b.x();
        core::Real sum_y = a.y() + b.y();
        core::Real sum_z = a.z() + b.z();
        
        return numeric::xyzVector< core::Real >( sum_x, sum_y, sum_z );
        
    }

    /// @brief    Difference from 2 Numeric XYZ Reals
    /// @details  Calculate the vector difference of 2 nymeric xyz vectors
    ///           utility function for the two methods below - maintains precision for core::Real
    numeric::xyzVector< core::Real >
    xyz_diff( numeric::xyzVector< core::Real > & a, numeric::xyzVector< core::Real > & b) {
        
        numeric::xyzVector< core::Real > diff;
        core::Real diff_x = a.x() - b.x();
        core::Real diff_y = a.y() - b.y();
        core::Real diff_z = a.z() - b.z();
        
        return numeric::xyzVector< core::Real >( diff_x, diff_y, diff_z );
        
    }
    
    /// @brief   Retrieve CA Coordiantes for a set of residues
    /// @details Given a set of residue numebrs, grab the CA xys coordinates at
    ///          that given residue position in the pose
    ///
    /// @return  Map for residue position to xyz coordinates
    std::map< core::Size, numeric::xyzVector< core::Real > >
    get_rsd_CAs( core::pose::Pose & pose, utility::vector1< core::Size > residues ) {
        
        using core::id::NamedAtomID;
        using numeric::xyzVector;
        
        // Create CA coords array and set dimension
        std::map< core::Size, numeric::xyzVector< core::Real > > CA_coords;
        
        // For each residue listed from the pose, grab coordinates, make
        // a new vector, and insert it into the iterable map
        for ( core::Size i = 1; i <= residues.size(); ++i ) {
            
            NamedAtomID id( "CA", residues[i] );
            xyzVector< core::Real > xyz( pose.xyz(id) );
            CA_coords.insert( std::pair< core::Real, numeric::xyzVector< core::Real > >( residues[i], xyz ));
        }
        
        // Return list of CA coordinates
        return CA_coords;
    }
    
    /// @brief   Retrieve CB Coordinates for a set of residues
    /// @details Given a set of residue numebrs, grab the CB xyz coordinates
    ///          at that given residue position
    ///
    /// @return  Map for residue position to xyz coords
    std::map< core::Size, numeric::xyzVector< core::Real > >
    get_rsd_CBs( core::pose::Pose & pose, utility::vector1< core::Size > residues ) {
        
        using core::id::NamedAtomID;
        using numeric::xyzVector;
        
        // Create CA coords array and set dimension
        std::map< core::Size, numeric::xyzVector< core::Real > > CB_coords;
        
        // For each residue listed from the pose, grab coordinates, make
        // a new vector, and insert it into the iterable map
        for ( core::Size i = 1; i <= residues.size(); ++i ) {
            
            NamedAtomID id( "CB", residues[i] );
            xyzVector< core::Real > xyz( pose.xyz(id) );
            CB_coords.insert( std::pair< core::Real, numeric::xyzVector< core::Real > >( residues[i], xyz ));
        }
        
        // Return list of CA coordinates
        return CB_coords;
    }
    
    /// @brief    Calculate net residue CA->COM vectors
    /// @details  Calculate the Normal vector as from the net residueCA-chain-COM vectors
    ///
    /// @param    pose
    ///             chainof interest
    /// @param    residueCOM
    ///             coordinates for CA of residue center of mass (for which CA coordinates will be grabbed)
    /// @param    residues
    ///             list of relevant residue positions
    numeric::xyzVector< core::Real >
    calc_net_CA_COM(
                    core::pose::Pose & pose,
                    numeric::xyzVector< core::Real > residueCOM,
                    utility::vector1< core::Size > residues
                    ) {
        
        using core::id::NamedAtomID;
        using numeric::xyzVector;
        
        // Grab residue CAs
        std::map< core::Size, numeric::xyzVector< core::Real > > residueCAs = get_rsd_CAs( pose, residues);
        
        // Initialize Sum
        numeric::xyzVector< core::Real > final_vec(0, 0, 0);
        
        // Grab each resCa_poseCOM vector and recursive add to the final net vector
        for ( core::Size i = 1; i <= residueCAs.size(); ++i ) {
            numeric::xyzVector< core::Real > xyz = xyz_diff( residueCOM, residueCAs.at( residues[i] ) );
            final_vec = xyz_sum( final_vec, xyz );
        }
        
        return final_vec;
    }
    
    /// @brief    Calculate net residue CA->CB Vectors
    /// @details  Calculate the Normal vector as net of CA->CB vectors
    ///
    /// @param    pose
    ///             pose of interest
    /// @param    residues
    ///             list of relevant residue positions
    numeric::xyzVector< core::Real >
    calc_net_CA_CB(
                   core::pose::Pose & pose,
                   utility::vector1< core::Size > residues
                   ) {
        
        using core::id::NamedAtomID;
        using numeric::xyzVector;
        
        // Grab residue CAs and CBs
        std::map< core::Size, numeric::xyzVector< core::Real > > residueCAs = get_rsd_CAs( pose, residues);
        std::map< core::Size, numeric::xyzVector< core::Real > > residueCBs = get_rsd_CBs( pose, residues);

        // Initialize Sum
        numeric::xyzVector< core::Real > final_vec(0, 0, 0);
        
        // Grab each resCa_poseCOM vector and recursive add to the final net vector
        for ( core::Size i = 1; i <= residueCAs.size(); ++i ) {
            numeric::xyzVector< core::Real > xyz = xyz_diff( residueCAs.at( residues [i] ), residueCBs.at( residues[i] ) );
            final_vec = xyz_sum( final_vec, xyz );
        }
        
        return final_vec;
        
    }
    
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
                   int const start,
                   int const stop
                   )
    {
        Vector center( 0.0 );
        for ( int i=start; i<=stop; ++i ) {
            if( !pose.residue( i ).is_protein()) {
                Vector ca_pos( pose.residue( i ).nbr_atom_xyz() );
                center += ca_pos;
            } else {
                Vector ca_pos( pose.residue( i ).atom( "CA" ).xyz() );
                center += ca_pos;
			}
        }
        center /= (stop-start+1);
        
        return center;
    }
    
    
    /// @brief      Residue Center of Mass
    /// @details    Calcualte the center of mass of a pose.
    ///
    /// @author     Monica Berrondo, Modified by Javier Castellanos and Rebecca Alford
    int
    residue_center_of_mass(
                           pose::Pose const & pose,
                           int const start,
                           int const stop
                           )
    {
        Vector center = center_of_mass(pose, start, stop );
        return return_nearest_residue( pose, start, stop, center );
    }
    
    /// @brief      Return nearest residue
    /// @details    Find the residue nearest some position passed in (normally a center of mass)
    ///
    /// @author     Monica Berrondo, Modified by Javier Castellanos and Rebecca Alford
    int
    return_nearest_residue(
                           pose::Pose const & pose,
                           int const begin,
                           int const end,
                           Vector center
                           )
    {
        Real min_dist = 9999.9;
        int res = 0;
        for ( int i=begin; i<=end; ++i )
        {
            Vector ca_pos;
            if( !pose.residue( i ).is_protein() ){
                ca_pos = pose.residue( i ).nbr_atom_xyz();
			} else {
                //Vector ca_pos( pose.residue( i ).atom( "CA" ).xyz() );
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
    
} // geometry
} // membrane
} // core

#endif // INCLUDED_core_membrane_geometry_util_cc

