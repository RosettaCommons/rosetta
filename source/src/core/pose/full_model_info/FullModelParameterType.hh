// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/full_model_info/FullModelParameterType.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_full_model_info_FullModelParameterType_HH
#define INCLUDED_core_pose_full_model_info_FullModelParameterType_HH

#include <core/pose/full_model_info/FullModelParameters.fwd.hh>

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
	// Note that this is also used to enum over chains and fixed_domain, which
	//  have the form [0,1,1,1,0,0,0,2,2,2,...].
	//
	// For these cases, we have an accessor function
	//    get_res_list_with_value( FIXED_DOMAIN, 3 ).
	//
	// -- rhiju, 2014.

	//
	// IF YOU ADD A TYPE, please also update initialize_parameters() in FullModelParameterType.cc
	//

	enum FullModelParameterType {
		NO_TYPE = 0,
		CALC_RMS,
		CHAIN, // Starts with 1 for first chain, then 2, ... note must agree with CUTPOINT_OPEN.
		CUTPOINT_OPEN, // Just points where a full pose would have fold-tree cutpoints.
		FIXED_DOMAIN,  // Non-zero over residues that are fixed (e.g., from input pdbs). Cannot delete or split these.
		EXTRA_MINIMIZE, // Specified by user for job -- can go into 'fixed domains'
		SAMPLE, // Whatever we can add/delete. Cannot include any FIXED_DOMAIN.

		// next-section: RNA-specific. Eventually expand.
		RNA_SYN_CHI,
		RNA_NORTH_SUGAR,
		RNA_SOUTH_SUGAR,
		RNA_TERMINAL, // base pairs on which new bases cannot stack/pair.

		// next-section: anything protein-specific

		// done.
		LAST_TYPE
	};

	void
	initialize_parameters( FullModelParameters & full_model_parameters );

} //full_model_info
} //pose
} //core

#endif
