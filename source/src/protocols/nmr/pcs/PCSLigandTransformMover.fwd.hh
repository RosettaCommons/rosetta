// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pcs/PCSLigandTransformMover.fwd.hh
/// @brief   forward declaration for PCSLigandTransformMover
/// @details last Modified: 06/02/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_nmr_pcs_PCSLigandTransformMover_FWD_HH
#define INCLUDED_protocols_nmr_pcs_PCSLigandTransformMover_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <platform/types.hh>

namespace protocols {
namespace nmr {
namespace pcs {

class AtomGridPoint;
class AtomGrid;
class PCSLigandTransformMover;

typedef utility::pointer::shared_ptr< PCSLigandTransformMover > PCSLigandTransformMoverOP;
typedef utility::pointer::shared_ptr< PCSLigandTransformMover const > PCSLigandTransformMoverCOP;
typedef utility::pointer::weak_ptr< PCSLigandTransformMover > PCSLigandTransformMoverAP;
typedef utility::pointer::weak_ptr< PCSLigandTransformMover const > PCSLigandTransformMoverCAP;

// Forward declaration of friend function
void pcs_lmmin(platform::Real const *par, int /*m_dat*/, void const *data, platform::Real *fvec, int */*info*/);

} // namespace pcs
} // namespace nmr
} // namespace protocols

#endif // INCLUDED_protocols_nmr_pcs_PCSLigandTransformMover_FWD_HH
