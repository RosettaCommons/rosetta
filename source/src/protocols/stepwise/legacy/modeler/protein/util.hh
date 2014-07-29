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


#ifndef INCLUDED_protocols_stepwise_legacy_modeling_protein_StepWiseProteinUtil_HH
#define INCLUDED_protocols_stepwise_legacy_modeling_protein_StepWiseProteinUtil_HH

#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/types.hh>

//Auto Headers
#include <utility/vector1.hh>

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace protein {

	utility::vector1< core::Size > const
	convert_to_working_res( utility::vector1< core::Size > const & res_vector,
													utility::vector1< core::Size > const & working_res ) ;

	void
	output_pose_list( utility::vector1< core::pose::PoseOP > pose_list,
										core::pose::PoseCOP native_pose,
										std::string const & silent_file,
										utility::vector1< core::Size > const & working_calc_rms_res	);

	void
	output_silent_struct( core::pose::Pose const & pose, core::pose::PoseCOP const & native_pose_op,
												std::string const & silent_file, std::string const & tag,
												core::io::silent::SilentFileDataOP sfd_in = 0);


	void
	output_silent_struct( core::pose::Pose const & pose, core::pose::PoseCOP const & native_pose_op,
												std::string const & silent_file, std::string const & tag,
												core::io::silent::SilentFileDataOP sfd_in,
												utility::vector1< core::Size > const & calc_rms_res	);

	core::Real get_pretend_psi_explicit( core::pose::Pose const & pose, core::Size const & res );

	core::Real get_pretend_phi_explicit( core::pose::Pose const & pose, core::Size const & res );

	bool
	is_close_chain_break( core::pose::Pose const & pose);

} //protein
} //modeler
} //legacy
} //stepwise
} //protocols


#endif
