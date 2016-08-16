// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/pose_inputters/PDBPoseInputter.fwd.hh
/// @brief  Forward declaration of the %PDBPoseInputter class for initializing Poses from .pdb or .pdb.gz files
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Steven Lewis (smlewi@gmail.com)


#ifndef INCLUDED_protocols_jd3_pose_inputters_PDBPoseInputter_FWD_HH
#define INCLUDED_protocols_jd3_pose_inputters_PDBPoseInputter_FWD_HH

//utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd3 {
namespace pose_inputters {

class PDBPoseInputter;

typedef utility::pointer::shared_ptr< PDBPoseInputter > PDBPoseInputterOP;
typedef utility::pointer::shared_ptr< PDBPoseInputter const > PDBPoseInputterCOP;


} // namespace pose_inputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_PDBPoseInputter_FWD_HH
