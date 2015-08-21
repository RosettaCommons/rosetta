// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PoseEvaluator
/// @brief PoseEvaluator
/// @details
///
///
/// @author Oliver Lange


// Unit Headers
#include <protocols/simple_filters/JumpEvaluator.hh>

// Package Headers

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>

//Auto Headers
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>


namespace protocols {
namespace simple_filters {

using namespace core;

/* It probably would be nice to change the interface away from using jump_nr towards
explicit naming of the Stubs one wants to evaluate

New interface would be :
NamedStubID up_stub
NamedStubID down_stub

that would enable us to get rid of the current distinction ( is_protein() ) which
switches between using the assumption that stub is N, CA, C or using the intrinsic atom-tree stubs

using the atom-tree stubs is a bit against the philosophy of this evaluator since its existence was actually
prompted by the problem that the atom-tree was using different atoms for the stub then expected and thus screwing
up the simulation. An Evaluator using the atom-tree stubs would not have detected this problem.

*/

JumpEvaluator::JumpEvaluator( pose::Pose const& native_pose, Size jump_nr ) :
	evaluation::SingleValuePoseEvaluator< core::Real >( "RT_"+ObjexxFCL::string_of( jump_nr ) )
	// jump_nr_( jump_nr )
{
	using namespace kinematics;
	kinematics::Edge jump_edge = native_pose.fold_tree().jump_edge( jump_nr );
	Size res1=jump_edge.start();
	Size res2=jump_edge.stop();

	// work out the stubID
	chemical::ResidueType const& rt1 ( native_pose.residue_type ( res1 ) );
	chemical::ResidueType const& rt2 ( native_pose.residue_type ( res2 ) );

	if ( jump_edge.has_atom_info() && !(rt1.is_protein() && rt2.is_protein()) ) { // if jump_atoms defined, get stubs directly from atom_tree. temporary fix to pass unit test, but check with Oliver!
		up_jump_atom_   = id::AtomID( rt1.atom_index( jump_edge.start_atom() ), res1 );
		down_jump_atom_ = id::AtomID( rt2.atom_index( jump_edge.stop_atom() ),  res2 );
		native_up_ = native_pose.conformation().atom_tree().atom( up_jump_atom_ ).get_stub();
		native_down_ = native_pose.conformation().atom_tree().atom( down_jump_atom_ ).get_stub();
	} else { // otherwise assume standard protein residue to residue jump
		id::AtomID a1( rt1.atom_index ("N") , res1 );
		id::AtomID a2( rt1.atom_index ("CA") , res1 );
		id::AtomID a3( rt1.atom_index ("C") , res1 );
		down_stub_ = id::StubID( a1, a2, a3 );

		id::AtomID b1( rt2.atom_index ("N") , res2 );
		id::AtomID b2( rt2.atom_index ("CA") , res2 );
		id::AtomID b3( rt2.atom_index ("C") , res2 );
		up_stub_ = id::StubID(b1, b2, b3 );

		native_up_ = native_pose.conformation().atom_tree().stub_from_id( up_stub_ );
		native_down_ = native_pose.conformation().atom_tree().stub_from_id( down_stub_ );
	}
}


core::Real
JumpEvaluator::apply(
	pose::Pose& pose
) const {

	kinematics::Stub up, down;
	if ( up_jump_atom_.valid() && down_jump_atom_.valid() ) { // jump_atoms defined
		up = pose.conformation().atom_tree().atom( up_jump_atom_ ).get_stub();
		down = pose.conformation().atom_tree().atom( down_jump_atom_ ).get_stub();
	} else {
		up = pose.conformation().atom_tree().stub_from_id( up_stub_ );
		down = pose.conformation().atom_tree().stub_from_id( down_stub_ );
	}
	kinematics::RT rt(up, down);

	kinematics::Stub test_down;
	rt.make_jump( native_up_, test_down );

	Real rms( 0.0 );
	for ( Size i=1; i<=3; i++ ) {
		Vector tv = test_down.build_fake_xyz( i );
		Vector nv = native_down_.build_fake_xyz( i );
		Vector d = nv-tv;
		rms += d.length();
	}
	return rms;
}

core::Size
JumpNrEvaluator::apply(
	pose::Pose& pose
) const {
	return pose.num_jump();
}


}
}
