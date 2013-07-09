// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file types.hh
/// @brief a couple enum types useful in SWA MonteCarlo applications
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_monte_carlo_types_hh
#define INCLUDED_protocols_swa_monte_carlo_types_hh

namespace protocols {
namespace swa {
namespace monte_carlo {

	enum MovingResidueCase { NO_CASE = 0, CHAIN_TERMINUS_5PRIME, CHAIN_TERMINUS_3PRIME, INTERNAL, FLOATING_BASE };

	enum AddOrDeleteChoice{ NO_ADD_OR_DELETE = 0, ADD, DELETE };

} // monte_carlo
} // swa
} // protocols

#endif
