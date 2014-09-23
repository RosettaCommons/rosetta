// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Environment.fwd.hh
/// @brief definition of the BrokeredEnvironment class
/// @author

#ifndef INCLUDED_protocols_environment_Environment_fwd_hh
#define INCLUDED_protocols_environment_Environment_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <boost/shared_ptr.hpp>

// Package headers

namespace protocols {
namespace environment {

class Environment;
typedef utility::pointer::shared_ptr< Environment > EnvironmentOP;
typedef utility::pointer::shared_ptr< Environment const > EnvironmentCOP;

typedef utility::pointer::weak_ptr< Environment > EnvironmentAP;
typedef utility::pointer::weak_ptr< Environment const > EnvironmentCAP;

} // environment
} // protocols

#endif //INCLUDED_protocols_moves_Environment_fwd_HH
