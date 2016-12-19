// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file core/chemical/PoseResidueTypeSet.fwd.hh
/// @brief A ResidueTypeSet which can be cached in the Pose
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_PoseResidueTypeSet_fwd_hh
#define INCLUDED_core_chemical_PoseResidueTypeSet_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace core {
namespace chemical {

class PoseResidueTypeSet;

typedef utility::pointer::shared_ptr< PoseResidueTypeSet > PoseResidueTypeSetOP;
typedef utility::pointer::shared_ptr< PoseResidueTypeSet const > PoseResidueTypeSetCOP;



} //core
} //chemical


#endif //INCLUDED_core_chemical_PoseResidueTypeSet_fwd_hh





