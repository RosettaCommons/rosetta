// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/abinitio/BrokerMain.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_abinitio_BrokerMain_hh
#define INCLUDED_protocols_abinitio_BrokerMain_hh

// *** IMPORTANT ***
// register_options_broker() must be called prior to any of the *_main() methods

namespace protocols {
namespace abinitio {

/// @brief Registers options that are relevant to the application. This function
/// must be called prior to Broker_main() or Loopbuild_Threading_main().
void register_options_broker();

/// @brief Initiates AbrelaxMover using the job distribution (jd2)
void Broker_main();

}
}

#endif  // INCLUDED_protocols_abinitio_BrokerMain_hh
