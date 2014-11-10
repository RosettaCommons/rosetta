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

#ifndef INCLUDED_protocols_membrane_metrics_hh
#define INCLUDED_protocols_membrane_metrics_hh

// Package Headers
#include <core/pose/Pose.fwd.hh> 
#include <core/types.hh> 

namespace protocols {
namespace membrane {

using namespace core::pose;
    
/// @brief Compute No Super Membrane RMSD
/// @details Compute the root mean squared deviation of the backbone atoms between
/// in the transmembrane spanning regions of a membrane protein
/// @throws If not a membrane protein
core::Real
mem_bb_rmsd_no_super( Pose & native_pose, Pose & pose );

/// @brief Compute No Super Membrane RMSD - all atom without super
/// @details Compute the root mean squared deviation of the backbone atoms between
/// in the transmembrane spanning regions of a membrane protein
/// @throws If not a membrane protein
core::Real
mem_all_atom_rmsd_no_super( Pose & native_pose, Pose & pose );
    
/// @brief Compute Membrane RMSD
/// @details Compute the root mean squared deviation of the backbone atoms between
/// in the transmembrane spanning regions of a membrane protein
core::Real
mem_bb_rmsd_with_super( Pose & native_pose, Pose & pose );
    
/// @brief Compute Membrane RMSD - all atom with superposiiton
/// @details Compute the root mean squared deviation of the backbone atoms between
/// in the transmembrane spanning regions of a membrane protein
core::Real
mem_all_atom_rmsd_with_super( Pose & native_pose, Pose & pose );
    
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_metrics_hh
