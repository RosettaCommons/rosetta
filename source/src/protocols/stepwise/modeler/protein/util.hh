// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_ProtocolUtil.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_protein_StepWiseProteinUtil_HH
#define INCLUDED_protocols_stepwise_protein_StepWiseProteinUtil_HH

#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/types.hh>

//Auto Headers
#include <utility/vector1.hh>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace protein {


	// This is used by create_alignment_id_map, but should be deprecated in favor of StepWisePoseAligner
 	void
 	setup_protein_backbone_atom_id_map( core::pose::Pose const & pose_1,
																			core::pose::Pose const & pose_2,
																			core::Size const base_res,
																			core::id::AtomID_Map< core::id::AtomID > & atom_ID_map);

// This is used by create_alignment_id_map, but should be deprecated in favor of StepWisePoseAligner
 	void
 	setup_protein_backbone_atom_id_map( core::pose::Pose const & pose_1,
																			core::pose::Pose const & pose_2,
																			core::Size const base_res,
																			core::Size const base_res2,
																			core::id::AtomID_Map< core::id::AtomID > & atom_ID_map);

	void
	figure_out_protein_modeling_info( core::pose::Pose const & pose,
																		core::Size const moving_res,
																		utility::vector1< core::Size > & moving_res_list );

	utility::vector1< core::Size >
	get_bridge_res( core::pose::Pose const & pose,
									utility::vector1< core::Size > const & moving_res_list /*working*/ );

	utility::vector1< core::Size >
	just_protein( utility::vector1< core::Size > const & res_list, core::pose::Pose const & pose );

	bool
	contains_protein( core::pose::Pose const & pose );

} //protein
} //modeler
} //stepwise
} //protocols


#endif
