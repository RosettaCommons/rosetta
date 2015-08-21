// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/HelixRotate.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/HelixRotate.hh>

// C/C++ headers
#include <iostream>
#include <string>

// Utility header
#include <basic/Tracer.hh>
#include <numeric/xyzVector.hh>

// Project headers
#include <core/conformation/Conformation.hh>
#include <core/id/NamedAtomID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/nonlocal/StarTreeBuilder.hh>

// Package headers
#include <protocols/moves/Mover.hh>

//Auto Headers
#include <utility/vector1.hh>
namespace protocols {
namespace nonlocal {

static thread_local basic::Tracer TR( "protocols.nonlocal.HelixRotate" );

HelixRotate::HelixRotate() {
	initialize(protocols::loops::Loop(), 0.0);
}

HelixRotate::HelixRotate(const protocols::loops::Loop& helix, double degrees) {
	initialize(helix, degrees);
}

void HelixRotate::initialize(const protocols::loops::Loop& helix, double degrees) {
	helix_ = helix;
	degrees_ = degrees;
}

void HelixRotate::apply(core::pose::Pose& pose) {
	using core::id::NamedAtomID;
	using core::kinematics::FoldTree;
	using core::kinematics::Jump;
	using core::kinematics::Stub;
	using numeric::xyzVector;
	using protocols::loops::Loops;
	using protocols::nonlocal::StarTreeBuilder;

	if ( !is_valid() ) {
		TR.Warning << "HelixRotate::apply() invoked with invalid or incomplete information." << std::endl;
		TR.Warning << "  helix_ => " << get_helix() << std::endl;
		TR.Warning << "  degrees_ => " << get_degrees() << std::endl;
		return;
	}

	// Retain a copy of the input fold tree, since we're responsible for restoring it
	FoldTree input_tree = pose.fold_tree();

	// Configure new kinematics
	Loops chunks;
	decompose_structure(pose.total_residue(), &chunks);

	StarTreeBuilder builder;
	builder.set_up(chunks, &pose);

	// Define the axis of translation
	xyzVector<double> axis, point;
	get_rotation_parameters(pose, &axis, &point);

	// Rotation about the axis
	unsigned jump_num = jump_containing_helix(chunks);
	Jump jump = pose.jump(jump_num);
	jump.rotation_by_axis(pose.conformation().upstream_jump_stub(jump_num), axis, point, get_degrees());
	pose.set_jump(jump_num, jump);

	// Restore input fold tree
	builder.tear_down(&pose);
	pose.fold_tree(input_tree);
}

void avg_ca_position(
	const core::pose::Pose& pose,
	const protocols::loops::Loop& region,
	numeric::xyzVector<double>* point
) {
	assert(point);

	point->zero();
	for ( unsigned i = region.start(); i <= region.stop(); ++i ) {
		(*point) += pose.xyz(core::id::NamedAtomID("CA", i));
	}

	(*point) /= region.length();
}

void HelixRotate::get_rotation_parameters(
	const core::pose::Pose& pose,
	numeric::xyzVector<double>* axis,
	numeric::xyzVector<double>* point
) const
{
	using core::id::NamedAtomID;
	using numeric::xyzVector;
	using protocols::loops::Loop;
	assert(axis);
	assert(point);

	// Define the point of rotation to be the average of the midpoint +/- 1 residue
	avg_ca_position(pose, Loop(helix_.midpoint() - 1, helix_.midpoint() + 1), point);

	if ( helix_.length() < 6 ) {
		*axis = pose.xyz(NamedAtomID("CA", helix_.stop())) - pose.xyz(NamedAtomID("CA", helix_.start()));
		return;
	}

	// Given a sufficient number of points, define the axis of rotation by the
	// average position of the first and last 3 CA atoms.
	xyzVector<double> a, b;
	avg_ca_position(pose, Loop(helix_.start(), helix_.start() + 2), &a);
	avg_ca_position(pose, Loop(helix_.stop() - 2, helix_.stop()), &b);
	*axis = b - a;
}

unsigned HelixRotate::jump_containing_helix(const protocols::loops::Loops& chunks) const {
	for ( unsigned i = 1; i <= chunks.num_loop(); ++i ) {
		if ( chunks[i].start() == helix_.start() ) {
			return i;
		}
	}
	return 0;  // invalid
}

void HelixRotate::decompose_structure(unsigned num_residues, protocols::loops::Loops* chunks) const {
	using protocols::loops::Loop;
	assert(chunks);
	assert(num_residues > 0);

	const unsigned start = get_helix().start();
	const unsigned stop = get_helix().stop();

	// Residues 1 to (helix - 1)
	if ( start > 1 ) {
		chunks->add_loop(Loop(1, start - 1));
	}

	// Helix
	chunks->add_loop(Loop(start, stop));

	// Residues (helix + 1) to end
	if ( stop < num_residues ) {
		chunks->add_loop(Loop(stop + 1, num_residues));
	}

	chunks->sequential_order();
}

bool HelixRotate::is_valid() const {
	return helix_.start() > 0 && helix_.start() < helix_.stop();
}

const protocols::loops::Loop& HelixRotate::get_helix() const {
	return helix_;
}

void HelixRotate::set_helix(const protocols::loops::Loop& helix) {
	helix_ = helix;
}

double HelixRotate::get_degrees() const {
	return degrees_;
}

void HelixRotate::set_degrees(double degrees) {
	degrees_ = degrees;
}

std::string HelixRotate::get_name() const {
	return "HelixRotate";
}

moves::MoverOP HelixRotate::fresh_instance() const {
	return moves::MoverOP( new HelixRotate() );
}

moves::MoverOP HelixRotate::clone() const {
	return moves::MoverOP( new HelixRotate(*this) );
}

}  // namespace nonlocal
}  // namespace protocols
