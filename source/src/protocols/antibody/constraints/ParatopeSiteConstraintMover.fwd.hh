// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody_design/ParatopeSiteConstraintMover.fwd.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_constraints_PARATOPESITECONSTRAINTMOVER_FWD_HH
#define INCLUDED_protocols_antibody_constraints_PARATOPESITECONSTRAINTMOVER_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

// C++ Headers
namespace protocols {
namespace antibody {
namespace constraints {


// Forward
class ParatopeSiteConstraintMover;

typedef utility::pointer::shared_ptr< ParatopeSiteConstraintMover > ParatopeSiteConstraintMoverOP;
typedef utility::pointer::shared_ptr< ParatopeSiteConstraintMover const > ParatopeSiteConstraintMoverCOP;


}
} //namespace antibody
} //namespace protocols

#endif //#ifndef INCLUDED_protocols/antibody_design_PARATOPESITECONSTRAINTMOVER_FWD_HH

