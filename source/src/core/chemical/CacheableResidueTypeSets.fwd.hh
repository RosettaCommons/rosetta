// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/CacheableResidueTypeSets.fwd.hh
/// @brief A (Pose-cacheable) container for ResidueTypeSets
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_CacheableResidueTypeSets_fwd_hh
#define INCLUDED_core_chemical_CacheableResidueTypeSets_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace core {
namespace chemical {

class CacheableResidueTypeSets;

typedef utility::pointer::shared_ptr< CacheableResidueTypeSets > CacheableResidueTypeSetsOP;
typedef utility::pointer::shared_ptr< CacheableResidueTypeSets const > CacheableResidueTypeSetsCOP;



} //core
} //chemical


#endif //INCLUDED_core_chemical_CacheableResidueTypeSets_fwd_hh





