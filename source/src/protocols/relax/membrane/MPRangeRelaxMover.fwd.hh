// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/relax/membrane/MPRangeRelaxMover.hh
/// @brief      Relaxes a membrane protein by relaxing in ranges
/// @details Relaxes a membrane protein by iteratively relaxing ranges of the protein;
///    No ramping required. Much faster than FastRelax and good for
///    large to very large proteins (tested up to 5250 residues)
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_relax_membrane_MPRangeRelaxMover_fwd_hh
#define INCLUDED_protocols_relax_membrane_MPRangeRelaxMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace relax {
namespace membrane {

class MPRangeRelaxMover;
typedef utility::pointer::shared_ptr< MPRangeRelaxMover > MPRangeRelaxMoverOP;
typedef utility::pointer::shared_ptr< MPRangeRelaxMover const > MPRangeRelaxMoverCOP;

} // membrane
} // relax
} // protocols

#endif // INCLUDED_protocols_relax_membrane_MPRangeRelaxMover_fwd_hh
