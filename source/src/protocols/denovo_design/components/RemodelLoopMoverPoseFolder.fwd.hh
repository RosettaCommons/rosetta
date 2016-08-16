// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/components/RemodelLoopMoverPoseFolder.fwd.hh
/// @brief Folds residues in a pose using RemodelLoopMover
/// @author Tom Linsky (tlinsky@uw.edu)


#ifndef INCLUDED_protocols_denovo_design_components_RemodelLoopMoverPoseFolder_fwd_hh
#define INCLUDED_protocols_denovo_design_components_RemodelLoopMoverPoseFolder_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace denovo_design {
namespace components {

class RemodelLoopMoverPoseFolder;

typedef utility::pointer::shared_ptr< RemodelLoopMoverPoseFolder > RemodelLoopMoverPoseFolderOP;
typedef utility::pointer::shared_ptr< RemodelLoopMoverPoseFolder const > RemodelLoopMoverPoseFolderCOP;



} //protocols
} //denovo_design
} //components


#endif //INCLUDED_protocols_denovo_design_components_RemodelLoopMoverPoseFolder_fwd_hh





