// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/membrane/metrics.hh
///
/// @brief   Metrics for membrane framework proteins
/// @details Metris for evaluating membrane protein models. Includes currently
///          a membrane RMSD metric w/ and w/o superimposition
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <protocols/membrane/metrics.hh> 

// Project Headers
#include <core/conformation/Conformation.hh> 
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh> 

#include <core/scoring/rms_util.hh> 
#include <core/scoring/rms_util.tmpl.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

namespace protocols {
namespace membrane {

using namespace core::pose;
using namespace core::scoring;
using namespace core::conformation::membrane;
    
/// @brief Compute No Super Membrane RMSD
/// @details Compute the root mean squared deviation of the backbone atoms between
/// in the transmembrane spanning regions of a membrane protein
/// @throws If not a membrane protein
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

/// @brief Compute No Super Membrane RMSD - all atom without super
/// @details Compute the root mean squared deviation of the backbone atoms between
/// in the transmembrane spanning regions of a membrane protein
/// @throws If not a membrane protein
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
        
/// @brief Compute Membrane RMSD
/// @details Compute the root mean squared deviation of the backbone atoms between
/// in the transmembrane spanning regions of a membrane protein
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

/// @brief Compute Membrane RMSD - all atom with superposiiton
/// @details Compute the root mean squared deviation of the backbone atoms between
/// in the transmembrane spanning regions of a membrane protein
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

} // membrane
} // protocols
