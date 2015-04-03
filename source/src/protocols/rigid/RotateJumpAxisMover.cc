// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rigid/RotateJumpAxisMover.cc
/// @brief RotateJumpAxisMover methods implemented
/// @author Steven Lewis

// Unit Headers
#include <protocols/rigid/RotateJumpAxisMover.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh> //get_anchor_and_root_atoms
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Jump.hh>

// Numeric Headers
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh> //degrees-radians

#include <numeric/random/random.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.rigid.RotateJumpAxisMover" );

namespace protocols {
namespace rigid {

void RotateJumpAxisMover::apply( core::pose::Pose & pose ){

	//first determine where the jump is
	core::kinematics::Edge const & rb_jump(pose.fold_tree().jump_edge(rb_jump_num_));
	core::Size const upstream_resid(rb_jump.start());
	core::Size const downstream_resid(rb_jump.stop());
	//these fail if for an underspecified jump (meaning, !edge.has_atom_info())
	//std::string const upstream_atom(pose.fold_tree().jump_edge(rb_jump_num_).upstream_atom());
	//std::string const downstream_atom(pose.fold_tree().jump_edge(rb_jump_num_).downstream_atom());
	core::conformation::Residue const & upstream_res(pose.residue(upstream_resid));
	core::conformation::Residue const & downstream_res(pose.residue(downstream_resid));

	//Use Phil's lookup to protect from !edge.has_atom_info()
	core::Size upstream_atomno, downstream_atomno;
	core::conformation::get_anchor_and_root_atoms(
		upstream_res,
		downstream_res,
		rb_jump,
		upstream_atomno,
		downstream_atomno);

	//calculate the rotation axis and angle
	//looking down the axis from the upstream to downstream atom, positive rotations are counterclockwise
	core::Vector axis( upstream_res.atom(upstream_atomno).xyz()//minus
										 - downstream_res.atom(downstream_atomno).xyz() );
	core::Angle angle(calc_angle());

	TR << angle << std::endl;
	numeric::xyzMatrix< core::Angle > rotation_matrix( numeric::rotation_matrix(axis, angle));

	//make new jump and apply to pose
	using core::kinematics::Stub;
	using core::kinematics::Jump;
	Stub const S1( pose.conformation().upstream_jump_stub( rb_jump_num_ ) );
	Stub const S2( pose.conformation().downstream_jump_stub( rb_jump_num_ ) );
	pose.set_jump(rb_jump_num_, Jump(S1, Stub(rotation_matrix * S2.M, S2.v )));
	return;
}//apply

std::string
RotateJumpAxisMover::get_name() const {
	return "RotateJumpAxisMover";
}

core::Angle RotateJumpAxisMover::calc_angle()
{	return numeric::conversions::radians(lower_angle_ + ((upper_angle_ - lower_angle_) * numeric::random::rg().uniform())); }

/// @details random angle constructor.  rb_jump_num is the number of the jump.  Magic numbers 180 and -179.9999999... maintain the uniform range.  I'm sure there's a better way to get [180, -180) but I can't figure out what it is.
RotateJumpAxisMover::RotateJumpAxisMover( core::Size const rb_jump_num )
	: moves::Mover(), rb_jump_num_(rb_jump_num), upper_angle_(180.0), lower_angle_(-179.9999999999999999999999999999999999999999999999)
{	Mover::type( "RotateJumpAxisMover" ); }

/// @details range of angles constructor - takes DEGREES not RADIANS.  rb_jump_num is the number of the jump.
RotateJumpAxisMover::RotateJumpAxisMover( core::Size const rb_jump_num, core::Angle const upper, core::Angle const lower )
	: moves::Mover(), rb_jump_num_(rb_jump_num), upper_angle_(upper), lower_angle_(lower)
{	Mover::type( "RotateJumpAxisMover" ); }

/// @details particular angle constructor - takes DEGREES not RADIANS.  rb_jump_num is the number of the jump.
RotateJumpAxisMover::RotateJumpAxisMover( core::Size const rb_jump_num, core::Angle const angle )
	: moves::Mover(), rb_jump_num_(rb_jump_num), upper_angle_(angle), lower_angle_(angle)
{	moves::Mover::type( "RotateJumpAxisMover" ); }

RotateJumpAxisMover::~RotateJumpAxisMover(){}

}//rigid
}//protocols
