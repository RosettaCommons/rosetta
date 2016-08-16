// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/abscript/StagePreparer.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/abinitio/abscript/StagePreparer.hh>

// Package headers
#include <core/environment/DofPassport.hh>
#include <protocols/environment/DofUnlock.hh>

// Project headers

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

static THREAD_LOCAL basic::Tracer tr( "protocols.environment.movers.StagePreparer", basic::t_info );

namespace protocols {
namespace abinitio {
namespace abscript {

using namespace core::environment;


} // abscript
} // abinitio
} // protocols
