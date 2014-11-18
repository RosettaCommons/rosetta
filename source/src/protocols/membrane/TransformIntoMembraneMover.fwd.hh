// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/TransformInfoMembraneMoverCreator.hh
/// @brief      Transform pose into membrane coordinates (Rosetta Scripts Hook)
/// @details	Requires a MembraneInfo object with all of its associated
///				information. This can be done by calling AddMembraneMover
///				beforehand. pose.conformation().is_membrane() should always
///				return true.
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_TransformIntoMembraneMover_fwd_hh
#define INCLUDED_protocols_membrane_TransformIntoMembraneMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {

class TransformIntoMembraneMover;
typedef utility::pointer::shared_ptr< TransformIntoMembraneMover > TransformIntoMembraneMoverOP;
typedef utility::pointer::shared_ptr< TransformIntoMembraneMover const > TransformIntoMembraneMoverCOP;

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_TransformIntoMembraneMover_fwd_hh
