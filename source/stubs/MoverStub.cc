// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/MoverStub.cc
/// @brief MoverStub methods implemented
/// @author


// Unit Headers
#include <protocols/moves/MoverStub.hh>


// Package Headers

// Project Headers
#include <core/pose/Pose.hh>



// Utility Headers
#include <core/util/Tracer.hh>
#include <core/types.hh>
#include <utility/exit.hh>
#include <numeric/random/random.hh>

// C++ Headers

using core::util::T;
using core::util::Error;
using core::util::Warning;

static core::util::Tracer TR( "protocols.moves.MoverStub" );
static numeric::random::RandomGenerator RG(MAKEMEANUMBERIFYOUNEEDME);

namespace protocols {
namespace moves {

/// @details
void MoverStub::apply( core::pose::Pose & pose ){

}//apply

/// @brief
MoverStub::MoverStub(
) : Mover()
{
	Mover::type( "MoverStub" );
}

MoverStub::~MoverStub(){}

}//moves
}//protocols

