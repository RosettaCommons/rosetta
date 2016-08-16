// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/denovo_design/movers/DisulfidizerMover.fwd.hh
/// @brief  DisulfidizerMover forward header
/// @author Tom Linsky (tlinsky@uw.edu) -- Adapting code from remodelmover into a mover
/// @author Gabe Rocklin (grocklin@uw.edu) -- Disulfide code


#ifndef INCLUDED_protocols_denovo_design_movers_DisulfidizerMover_fwd_hh
#define INCLUDED_protocols_denovo_design_movers_DisulfidizerMover_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace denovo_design {

// Forward
class DisulfidizerMover;

// Types
typedef  utility::pointer::shared_ptr< DisulfidizerMover >  DisulfidizerMoverOP;
typedef  utility::pointer::shared_ptr< DisulfidizerMover const >  DisulfidizerMoverCOP;

typedef  utility::pointer::weak_ptr< DisulfidizerMover >  DisulfidizerMoverAP;
typedef  utility::pointer::weak_ptr< DisulfidizerMover const >  DisulfidizerMoverCAP;


} // namespace denovo_design
} // namespace protocols

#endif
