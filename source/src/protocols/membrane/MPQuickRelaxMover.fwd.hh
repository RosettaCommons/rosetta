// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief      Does a quick relax of a membrane protein
/// @details Uses SmallMover and ShearMover with adjustable maximum dihedral
///    angle changes, then repacking and a single round of minimization
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_MPQuickRelaxMover_fwd_hh
#define INCLUDED_protocols_membrane_MPQuickRelaxMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {

class MPQuickRelaxMover;
typedef utility::pointer::shared_ptr< MPQuickRelaxMover > MPQuickRelaxMoverOP;
typedef utility::pointer::shared_ptr< MPQuickRelaxMover const > MPQuickRelaxMoverCOP;

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_MPQuickRelaxMover_fwd_hh
