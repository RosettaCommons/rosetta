// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/full_model_info/FullModelParameterType.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_full_model_info_FullModelParameterType_HH
#define INCLUDED_core_pose_full_model_info_FullModelParameterType_HH

#include <core/pose/full_model_info/FullModelParameters.fwd.hh>
#include <string>

namespace core {
namespace pose {
namespace full_model_info {

// During StepWiseMonteCarlo, we need the pose to carry information
//  on what subset of the full model it has modeled.
//
// This FullModelParameterType is meant to be a key into a map of useful
//  residue number lists, all numbered so that the eventual full-length model
//  is 1-N (i.e., Rosetta numbering for some full-length pose).
//
// Example:  get_res_list( CUTPOINT_OPEN )    --> [5,10].
//           value_at_res( 5, CUTPOINT_OPEN ) --> 1
//           value_at_res( 6, CUTPOINT_OPEN ) --> 0
//
// (The actual variable is stored as [0,0,0,0,1,0,0,0,0,1,..] as well as [5,10]. )
//
// One constraint -- these really should be 'permanent' -- parameters that won't
//  change during a run. For example, the 'bridge_res' for cutpoint closure is
//  not appropriate to store here, as those values change from move to move.
//
// Note that this is also used to enum over fixed_domain, input_domain, etc., which
//  have the form [0,1,1,1,0,0,0,2,2,2,...].
//
// For these cases, we have an accessor function
//    get_res_list_with_value( FIXED_DOMAIN, 3 ).
//
// -- rhiju, 2014.


// IF YOU ADD A TYPE, please also update initialize_parameters() in FullModelParameterType.cc,
//  as well as  initialize_full_model_parameter_type_name().

enum FullModelParameterType {
	NO_TYPE = 0,
	CALC_RMS,
	CUTPOINT_OPEN,  // Just points where a full pose would have fold-tree cutpoints.
	DOCK_DOMAIN,    // 0,1,2,... for separate domains that can be docked/undocked.
	EXTRA_MINIMIZE, // Specified by user for job -- can carve out of fixed domains
	EXTRA_MINIMIZE_JUMP, // Pairs of residues connected by jumps that should be minimized
	FIXED_DOMAIN,   // Non-zero over residues that are fixed. Can move/minimize relative to each other.
	INPUT_DOMAIN,   // from input PDBs, typically a superset of FIXED_DOMAIN
	JUMP,           // User-specified jump pairs, e.g., for fold-tree setup or docking contacts.
	PREFERRED_ROOT, // good positions for roots.
	SAMPLE,         // Whatever we can add/delete. Cannot include any FIXED_DOMAIN.
	WORKING,        // Everything that might be modeled [SAMPLE + FIXED_DOMAIN]
	// next-section: RNA-specific. Eventually expand.
	RNA_BULGE, // nucleotides that should not be instantiated.
	ALIGNMENT_ANCHOR_RES, // residue(s) that are part of the fixed domain to use for alignment
	RNA_SAMPLE_SUGAR, // residues that might be in fixed domains but allow sugar to move.
	RNA_SYN_CHI,
	RNA_ANTI_CHI,
	RNA_NORTH_SUGAR,
	RNA_SOUTH_SUGAR,
	RNA_TERMINAL, // base pairs on which new bases cannot stack
	RNA_BLOCK_STACK_ABOVE, // base pairs on which new bases cannot stack 'above' (3' direction)
	RNA_BLOCK_STACK_BELOW, // base pairs on which new bases cannot stack 'below' (5' direction)
	// next-section: anything protein-specific
	DISULFIDE,

	// done.
	LAST_TYPE
};

void
initialize_parameters( FullModelParameters & full_model_parameters );

std::string to_string( FullModelParameterType const & type );

FullModelParameterType
full_model_parameter_type_from_string( std::string const & name );

} //full_model_info
} //pose
} //core

#endif
