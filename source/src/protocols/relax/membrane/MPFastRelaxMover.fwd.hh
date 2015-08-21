// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/relax/membrane/MPFastRelaxMover.fwd.hh
///
/// @brief      Membrane Fast Relax Protocol - Relax with minimization of mem position
/// @details Apply the standard fast relax protocol. Enable minimization of the memrbane
///             jump and relax from the center of mass. Also use the smoothed
///             full atom membrane energy function.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified: 12/2/14

#ifndef INCLUDED_protocols_relax_membrane_MPFastRelaxMover_fwd_hh
#define INCLUDED_protocols_relax_membrane_MPFastRelaxMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace relax {
namespace membrane {

class MPFastRelaxMover;
typedef utility::pointer::shared_ptr< MPFastRelaxMover > MPFastRelaxMoverOP;
typedef utility::pointer::shared_ptr< MPFastRelaxMover const > MPFastRelaxMoverCOP;

} // membrane
} // relax
} // protocols

#endif // INCLUDED_protocols_relax_membrane_MPFastRelaxMover_fwd_hh
