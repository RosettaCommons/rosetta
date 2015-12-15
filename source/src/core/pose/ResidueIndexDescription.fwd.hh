// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/pose/ResidueIndexDescription.fwd.hh
/// @brief  Two classes designed to hold data neceassary to describe a residue in a Pose,
///         which may come from a text file, e.g., and to resolve that data into an actual
///         residue index when a Pose becomes available (which is likely not at the time
///         that the file is read) and to throw an exception if the index cannot be resolved.
/// @author Brian D. Weitzner
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_pose_ResidueIndexDescription_FWD_HH
#define INCLUDED_core_pose_ResidueIndexDescription_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pose {

class ResidueIndexDescription;
class ResidueIndexDescriptionFromFile;

typedef utility::pointer::shared_ptr< ResidueIndexDescription > ResidueIndexDescriptionOP;
typedef utility::pointer::shared_ptr< ResidueIndexDescription const > ResidueIndexDescriptionCOP;

typedef utility::pointer::shared_ptr< ResidueIndexDescriptionFromFile > ResidueIndexDescriptionFromFileOP;
typedef utility::pointer::shared_ptr< ResidueIndexDescriptionFromFile const > ResidueIndexDescriptionFromFileCOP;

}
}

#endif
