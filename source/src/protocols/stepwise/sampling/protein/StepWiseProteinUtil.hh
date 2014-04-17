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
#include <core/fragment/ConstantLengthFragSet.fwd.hh>

//Auto Headers
#include <utility/vector1.hh>

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

	core::Real
	get_rotamer_angle( core::Size const & i, core::Size const & N_SAMPLE );

	void
	output_silent_struct( core::pose::Pose const & pose, core::pose::PoseCOP const & native_pose_op,
												std::string const & silent_file, std::string const & tag,
												core::io::silent::SilentFileDataOP sfd_in = 0);


	void
	output_silent_struct( core::pose::Pose const & pose, core::pose::PoseCOP const & native_pose_op,
												std::string const & silent_file, std::string const & tag,
												core::io::silent::SilentFileDataOP sfd_in,
												utility::vector1< core::Size > const & calc_rms_res	);

	void
	remove_end_variants( core::pose::Pose & pose );

	core::Real get_pretend_psi_explicit( core::pose::Pose const & pose, core::Size const & res );

	core::Real get_pretend_phi_explicit( core::pose::Pose const & pose, core::Size const & res );

 	void
 	setup_protein_CA_atom_id_map( core::pose::Pose const & pose_1,
															  core::pose::Pose const & pose_2,
																core::Size const base_res,
																core::id::AtomID_Map< core::id::AtomID > & atom_ID_map);


 	void
 	setup_protein_backbone_atom_id_map( core::pose::Pose const & pose_1,
																			core::pose::Pose const & pose_2,
																			core::Size const base_res,
																			core::id::AtomID_Map< core::id::AtomID > & atom_ID_map);

 	void
 	setup_protein_backbone_atom_id_map( core::pose::Pose const & pose_1,
																			core::pose::Pose const & pose_2,
																			core::Size const base_res,
																			core::Size const base_res2,
																			core::id::AtomID_Map< core::id::AtomID > & atom_ID_map);

} //protein
} //sampling
} //stepwise
} //protocols


#endif
