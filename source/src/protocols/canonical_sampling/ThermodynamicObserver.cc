// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/canonical_sampling/ThermodynamicObserver.cc
/// @brief ThermodynamicObserver methods implemented
/// @author


// Unit Headers
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>


// Package Headers

// Project Headers


// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>

#include <utility/vector1.hh>


// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;

//static basic::Tracer TR( "protocols.canonical_sampling.ThermodynamicObserver" );

namespace protocols {
namespace canonical_sampling {

/// @brief
ThermodynamicObserver::ThermodynamicObserver(
) : Mover()
{
	Mover::type( "ThermodynamicObserver" );
}

ThermodynamicObserver::~ThermodynamicObserver() {}

} //moves
} //protocols

