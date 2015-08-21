// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/SpinAroundPartnerMover.fwd.hh
/// @brief      Spins the downstream partner around the upstream partner
/// @details Spins the downstream partner around the upstream partner in the
///    membrane to probe all kinds of interfaces. Both embedding normals
///    are approximately conserved, i.e. the partners aren't flipped
///    in the membrane.
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_SpinAroundPartnerMover_fwd_hh
#define INCLUDED_protocols_membrane_SpinAroundPartnerMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {

class SpinAroundPartnerMover;
typedef utility::pointer::shared_ptr< SpinAroundPartnerMover > SpinAroundPartnerMoverOP;
typedef utility::pointer::shared_ptr< SpinAroundPartnerMover const > SpinAroundPartnerMoverCOP;

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_SpinAroundPartnerMover_fwd_hh
