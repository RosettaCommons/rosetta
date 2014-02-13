// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file EnvClaimBroker.fwd.hh
/// @brief definition of the EnvClaimBroker class
/// @author

#ifndef INCLUDED_protocols_environment_EnvClaimBroker_fwd_hh
#define INCLUDED_protocols_environment_EnvClaimBroker_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <boost/shared_ptr.hpp>

// Package headers

namespace protocols {
namespace environment {

class EnvClaimBroker;
typedef utility::pointer::owning_ptr< EnvClaimBroker > EnvClaimBrokerOP;
typedef utility::pointer::owning_ptr< EnvClaimBroker const > EnvClaimBrokerCOP;

typedef utility::pointer::access_ptr< EnvClaimBroker > EnvClaimBrokerAP;
typedef utility::pointer::access_ptr< EnvClaimBroker const > EnvClaimBrokerCAP;

} // environment
} // protocols

#endif //INCLUDED_protocols_moves_EnvClaimBroker_fwd_hh
