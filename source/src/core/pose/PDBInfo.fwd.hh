// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/PDBInfo.fwd.hh
/// @brief  fwd declaration of classes defined in pose/PDBInfo
/// @author Steven Lewis (smlewi@gmail.com)


#ifndef INCLUDED_core_pose_PDBInfo_fwd_hh
#define INCLUDED_core_pose_PDBInfo_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.hh>

// C++ Headers
#include <string>

namespace core {
namespace pose {

class PDBInfo;

typedef utility::pointer::shared_ptr< PDBInfo > PDBInfoOP;
typedef utility::pointer::shared_ptr< PDBInfo const > PDBInfoCOP;

typedef  std::pair< std::string, std::string > ChainSegID;

} // namespace pose
} // namespace core

#endif
