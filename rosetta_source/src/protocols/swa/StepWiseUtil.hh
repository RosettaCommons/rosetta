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


#ifndef INCLUDED_protocols_swa_StepWiseUtil_hh
#define INCLUDED_protocols_swa_StepWiseUtil_hh

#include <core/pose/Pose.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh> // needed for default initialization of SilentFileDataOP
#include <core/types.hh>
#include <string>

namespace protocols {
namespace swa {

core::Real
get_rotamer_angle( core::Size const & i, core::Size const & N_SAMPLE );

void
output_silent_struct( core::pose::Pose const & pose, core::pose::PoseCOP const & native_pose_op,
											std::string const & silent_file, std::string const & tag,
											core::io::silent::SilentFileDataOP sfd_in = 0);

void
remove_end_variants( core::pose::Pose & pose );


core::Real get_pretend_psi_explicit( core::pose::Pose const & pose, core::Size const & res );

core::Real get_pretend_phi_explicit( core::pose::Pose const & pose, core::Size const & res );


} //swa
} // protocols

#endif
