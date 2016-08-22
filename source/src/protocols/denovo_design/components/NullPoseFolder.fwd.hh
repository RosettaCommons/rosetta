// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/components/NullPoseFolder.fwd.hh
/// @brief Pose folder that does nothing
/// @author Tom Linsky (tlinsky@uw.edu)


#ifndef INCLUDED_protocols_denovo_design_components_NullPoseFolder_fwd_hh
#define INCLUDED_protocols_denovo_design_components_NullPoseFolder_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace denovo_design {
namespace components {

class NullPoseFolder;

typedef utility::pointer::shared_ptr< NullPoseFolder > NullPoseFolderOP;
typedef utility::pointer::shared_ptr< NullPoseFolder const > NullPoseFolderCOP;



} //protocols
} //denovo_design
} //components


#endif //INCLUDED_protocols_denovo_design_components_NullPoseFolder_fwd_hh





