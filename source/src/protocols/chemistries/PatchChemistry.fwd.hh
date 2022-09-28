// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/chemistries/PatchChemistry.fwd.hh
/// @brief  Forward declaration of a class that applies Rosetta patches
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_chemistries_PatchChemistry_FWD_HH
#define INCLUDED_protocols_chemistries_PatchChemistry_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace chemistries {

class PatchChemistry;

typedef utility::pointer::shared_ptr< PatchChemistry > PatchChemistryOP;
typedef utility::pointer::shared_ptr< PatchChemistry const > PatchChemistryCOP;

} //namespace chemistries
} //namespace protocols

#endif
