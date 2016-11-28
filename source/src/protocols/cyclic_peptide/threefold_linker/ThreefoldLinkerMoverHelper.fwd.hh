// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/threefold_linker/ThreefoldLinkerMoverHelper.fwd.hh
/// @brief A base class for helper objects that the ThreefoldLinkerMover uses to set up specific types
/// of threefold linkers.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


#ifndef INCLUDED_protocols_cyclic_peptide_threefold_linker_ThreefoldLinkerMoverHelper_fwd_hh
#define INCLUDED_protocols_cyclic_peptide_threefold_linker_ThreefoldLinkerMoverHelper_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace protocols {
namespace cyclic_peptide {
namespace threefold_linker {

class ThreefoldLinkerMoverHelper;

typedef utility::pointer::shared_ptr< ThreefoldLinkerMoverHelper > ThreefoldLinkerMoverHelperOP;
typedef utility::pointer::shared_ptr< ThreefoldLinkerMoverHelper const > ThreefoldLinkerMoverHelperCOP;

} //threefold_linker
} //protocols
} //cyclic_peptide

#endif //INCLUDED_protocols_cyclic_peptide_threefold_linker_ThreefoldLinkerMoverHelper_fwd_hh
