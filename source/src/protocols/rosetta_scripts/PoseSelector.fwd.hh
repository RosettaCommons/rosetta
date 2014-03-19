// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/rosetta_scripts/PoseSelector.fwd.hh
/// @brief  Forward declarations for PoseSelector
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifndef INCLUDED_protocols_rosetta_scripts_PoseSelector_fwd_hh
#define INCLUDED_protocols_rosetta_scripts_PoseSelector_fwd_hh

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace rosetta_scripts {

class PoseSelector;
typedef utility::pointer::owning_ptr< PoseSelector > PoseSelectorOP;
typedef utility::pointer::owning_ptr< PoseSelector const > PoseSelectorCOP;

} // rosetta_scripts
} // protocols

#endif //INCLUDED_protocols_rosetta_scripts_PoseSelector_fwd_hh
