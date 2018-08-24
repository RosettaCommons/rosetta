// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/pose/ResidueIndexDescription.fwd.hh
/// @brief  Classes designed to hold data neceassary to describe a residue in a Pose,
///         which may come from a text file, e.g., and to resolve that data into an actual
///         residue index when a Pose becomes available (which is likely not at the time
///         that the file is read) and to throw an exception if the index cannot be resolved.
/// @author Brian D. Weitzner
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_pose_ResidueIndexDescription_FWD_HH
#define INCLUDED_core_pose_ResidueIndexDescription_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pose {

class RID_Source;
class ResidueIndexDescription;

typedef utility::pointer::shared_ptr< RID_Source > RID_SourceOP;
typedef utility::pointer::shared_ptr< RID_Source const > RID_SourceCOP;

typedef utility::pointer::shared_ptr< ResidueIndexDescription > ResidueIndexDescriptionOP;
typedef utility::pointer::shared_ptr< ResidueIndexDescription const > ResidueIndexDescriptionCOP;

}
}

#endif
