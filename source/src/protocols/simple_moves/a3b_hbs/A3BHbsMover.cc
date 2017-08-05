// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/a3b_hbs/A3BHbsMover.cc
/// @brief A3BHbsMover methods implemented
/// @author Andrew Watkins, amw579@nyu.edu

// Unit Headers
#include <protocols/simple_moves/a3b_hbs/A3BHbsMover.hh>
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
// Utility Headers
#include <numeric/xyz.functions.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
#include <core/types.hh>

// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.a3b_hbs.A3BHbsMover" );


using namespace core;
using namespace conformation;
using namespace chemical;
using namespace core::id;

namespace protocols {
namespace simple_moves {
namespace a3b_hbs {
// awatkins: This should NEVER be called on residues that actually have HBS atom labels, fixing... 7/14/13
/*
kdrew: the helper function moves a single hbs with the given move
pose: pose to make change to
hbs_pre_pos: position of first residue of hbs to make change to
phi_angle: value to change phi angle to
psi_angle: value to change psi angle to
*/
void A3BHbsMover::apply( core::pose::Pose & pose )
{
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	TR << "current pos ("<< seq_pos_<< ") PHI: " << pose.phi(seq_pos_) << std::endl;
	pose.set_phi(seq_pos_, phi_angle_);
	TR << "new pos ("<< seq_pos_<< ") PHI: " << pose.phi(seq_pos_) << std::endl;
	//pose.dump_pdb( "rosetta_out_hbs_pre.pdb" );
	TR << "current pos ("<< seq_pos_<< ") PSI: " << pose.psi(seq_pos_) << std::endl;
	pose.set_psi(seq_pos_, psi_angle_);
	TR << "new pos ("<< seq_pos_<< ") PSI: " << pose.psi(seq_pos_) << std::endl;
	//pose.dump_pdb( "rosetta_out_hbs_pre_pos_t.pdb" );

}

std::string
A3BHbsMover::get_name() const {
	return "HbsMover";
}

/// @brief
A3BHbsMover::A3BHbsMover(
	core::Size seq_position
): Mover(), seq_pos_(seq_position)
{
	Mover::type( "A3BHbsMover" );
}

A3BHbsMover::A3BHbsMover(
	core::Size seq_position,
	core::Real phi_angle,
	core::Real psi_angle
): Mover(), seq_pos_(seq_position), phi_angle_(phi_angle), psi_angle_(psi_angle)
{
	Mover::type( "A3BHbsMover" );

}

A3BHbsMover::~A3BHbsMover(){}

}//hbs
}//simple_moves
}//protocols

