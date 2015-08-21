// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/membrane/TransformIntoMembraneMover.fwd.hh
/// @brief  Transform a pose into a membrane coordinate frame
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// @author     JKLeman (julia.koehler1982@gmail.com)
/// Last Modified: 6/11/15
/// #RosettaMPMover

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
