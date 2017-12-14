// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/oop/OopPuckMover.cc
/// @brief OopPuckMover methods implemented
/// @author Kevin Drew, kdrew@nyu.edu

// Unit Headers
#include <protocols/simple_moves/oop/OopPuckMover.fwd.hh>
#include <protocols/simple_moves/oop/OopPuckMover.hh>
#include <protocols/simple_moves/oop/OopMover.hh>
#include <protocols/simple_moves/oop/OopPatcher.hh>
// Package Headers

// Project Headers
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
// Random number generator
#include <numeric/random/random.hh>
// Utility Headers
#include <numeric/xyz.functions.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
#include <core/types.hh>

// C++ Headers

using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.simple_moves.oop.OopPuckMover" );


using namespace core;
using namespace conformation;
using namespace chemical;
using namespace core::id;

//kdrew: defining constants
static const Real OOP_PUCK_PLUS_PHI = -131.27;
static const Real OOP_PUCK_PLUS_PSI = -10.38;

static const Real OOP_PUCK_MINUS_PHI = -147.05;
static const Real OOP_PUCK_MINUS_PSI = -36.90;


static const Real OOP_D_PUCK_PLUS_PHI = 146.05;
static const Real OOP_D_PUCK_PLUS_PSI = 7.40;

static const Real OOP_D_PUCK_MINUS_PHI = 151.90;
static const Real OOP_D_PUCK_MINUS_PSI = 30.40;


namespace protocols {
namespace simple_moves {
namespace oop {


OopPuckPlusMover::OopPuckPlusMover( core::Size oop_seq_position ): OopMover( oop_seq_position, OOP_PUCK_PLUS_PHI, OOP_PUCK_PLUS_PSI )
{
	OopMover::type( "OopPuckPlusMover" );
}

OopPuckPlusMover::~OopPuckPlusMover()= default;

std::string
OopPuckPlusMover::get_name() const {
	return "OopPuckPlusMover";
}

OopPuckMinusMover::OopPuckMinusMover( core::Size oop_seq_position ): OopMover( oop_seq_position, OOP_PUCK_MINUS_PHI, OOP_PUCK_MINUS_PSI )
{
	OopMover::type( "OopPuckMinusMover" );
}

OopPuckMinusMover::~OopPuckMinusMover()= default;

std::string
OopPuckMinusMover::get_name() const {
	return "OopPuckMinusMover";
}

OopDPuckPlusMover::OopDPuckPlusMover( core::Size oop_seq_position ): OopMover( oop_seq_position, OOP_D_PUCK_PLUS_PHI, OOP_D_PUCK_PLUS_PSI )
{
	OopMover::type( "OopDPuckPlusMover" );
}

OopDPuckPlusMover::~OopDPuckPlusMover()= default;

std::string
OopDPuckPlusMover::get_name() const {
	return "OopDPuckPlusMover";
}

OopDPuckMinusMover::OopDPuckMinusMover( core::Size oop_seq_position ): OopMover( oop_seq_position, OOP_D_PUCK_MINUS_PHI, OOP_D_PUCK_MINUS_PSI )
{
	OopMover::type( "OopDPuckMinusMover" );
}

OopDPuckMinusMover::~OopDPuckMinusMover()= default;

std::string
OopDPuckMinusMover::get_name() const {
	return "OopDPuckMinusMover";
}


}//oop
}//simple_moves
}//protocols

