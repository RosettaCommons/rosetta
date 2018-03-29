// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/monte_carlo/MonteCarloInterface.fwd.hh
/// @brief A MonteCarlo object for optimizing the interface dG as defined using InterfaceAnalyzer. The dG and Total energy can be weighted.  This is so that the interface energy itself can be optimized through a protocol.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_monte_carlo_MonteCarloInterface_fwd_hh
#define INCLUDED_protocols_monte_carlo_MonteCarloInterface_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace monte_carlo {

class MonteCarloInterface;

typedef utility::pointer::shared_ptr< MonteCarloInterface > MonteCarloInterfaceOP;
typedef utility::pointer::shared_ptr< MonteCarloInterface const > MonteCarloInterfaceCOP;

} //protocols
} //monte_carlo

#endif //INCLUDED_protocols_monte_carlo_MonteCarloInterface_fwd_hh
