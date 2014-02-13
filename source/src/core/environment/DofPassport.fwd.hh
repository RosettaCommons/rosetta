// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DofPassport.fwd.hh
/// @brief definition of the DofPassport class
/// @author

#ifndef INCLUDED_core_environment_DofPassport_fwd_hh
#define INCLUDED_core_environment_DofPassport_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <boost/shared_ptr.hpp>

// Package headers

namespace core {
namespace environment {

class DofPassport;
typedef utility::pointer::owning_ptr< DofPassport > DofPassportOP;
typedef utility::pointer::owning_ptr< DofPassport const > DofPassportCOP;

} // core
} // environment

#endif //INCLUDED_core_moves_DofPassport_fwd_HH
