// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rigid/RigidBodyMotionMover.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/rigid/RigidBodyMotionMover.hh>

// C/C++ headers
#include <cmath>
#include <iostream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/rigid.OptionKeys.gen.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>

// Project headers
#include <core/id/NamedAtomID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Package headers
#include <protocols/moves/Mover.hh>

//Auto Headers
#include <utility/vector1.hh>
#ifdef WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#endif


namespace protocols {
namespace rigid {

static THREAD_LOCAL basic::Tracer TR( "protocols.moves.RigidBodyMotionMover" );

double angle_between(const numeric::xyzVector<double>& a, const numeric::xyzVector<double>& b) {
	double radians = std::acos(a.dot(b) / (a.length() * b.length()));
	return radians * 180 / M_PI;
}

RigidBodyMotionMover::RigidBodyMotionMover(const core::kinematics::FoldTree& tree) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// update chunks
	set_fold_tree(tree);

	// retrieve defaults from options system
	set_magnitude_rotation(option[OptionKeys::rigid::rotation]());
	set_magnitude_translation(option[OptionKeys::rigid::translation]());
	set_chainbreak_bias(option[OptionKeys::rigid::chainbreak_bias]());
}

void RigidBodyMotionMover::apply(core::pose::Pose& pose) {
	using core::kinematics::Jump;
	using numeric::xyzVector;
	using std::endl;

	if ( chunks_.size() < 2 ) {
		TR << "No sensible action possible-- num_chunks=" << chunks_.size() << endl;
		return;
	} else if ( fold_tree() != pose.fold_tree() ) {
		TR << "Mismatching fold trees" << endl;
		TR << "Input: " << fold_tree() << endl;
		TR << "Apply: " << pose.fold_tree() << endl;
		return;
	}

	Size i = numeric::random::random_range(1, chunks_.size());
	Jump jump = pose.jump(i);

	// rotation
	jump.set_rb_delta(Jump::ROT_X, 1, numeric::random::gaussian() * magnitude_rotation());
	jump.set_rb_delta(Jump::ROT_Y, 1, numeric::random::gaussian() * magnitude_rotation());
	jump.set_rb_delta(Jump::ROT_Z, 1, numeric::random::gaussian() * magnitude_rotation());

	// translation
	xyzVector<double> bias;
	compute_bias(i, pose, &bias);

	double tx = numeric::random::gaussian() * magnitude_translation() + bias.x();
	double ty = numeric::random::gaussian() * magnitude_translation() + bias.y();
	double tz = numeric::random::gaussian() * magnitude_translation() + bias.z();

	jump.set_rb_delta(Jump::TRANS_X, 1, tx);
	jump.set_rb_delta(Jump::TRANS_Y, 1, ty);
	jump.set_rb_delta(Jump::TRANS_Z, 1, tz);

	// update the jump
	jump.fold_in_rb_deltas();
	pose.set_jump(i, jump);
}

void RigidBodyMotionMover::compute_bias(Size i, const core::pose::Pose& pose, numeric::xyzVector<double>* bias) const {
	using core::id::NamedAtomID;
	using numeric::xyzVector;
	using protocols::loops::Loop;
	debug_assert(bias);

	bool has_prev = i > 1;
	bool has_next = i < chunks_.size();

	const Loop& chunk = chunks_[i];

	if ( has_prev && has_next ) {
		const Loop& prev = chunks_[i - 1];
		const Loop& next = chunks_[i + 1];

		xyzVector<double> a = pose.xyz(NamedAtomID("CA", chunk.start())) - pose.xyz(NamedAtomID("CA", prev.stop()));
		xyzVector<double> b = pose.xyz(NamedAtomID("CA", chunk.stop())) - pose.xyz(NamedAtomID("CA", next.start()));

		// If the absolute value of the angle between vectors is too great,
		// we're better served biasing translation toward one of the endpoints
		double degrees = std::abs(angle_between(a, b));
		if ( degrees < 90 ) {
			*bias = (a + b) / 2.0;
		} else {
			*bias = (numeric::random::uniform() < 0.9) ? a : b;
		}
	} else if ( has_prev ) {
		const Loop& prev = chunks_[i - 1];
		*bias = pose.xyz(NamedAtomID("CA", chunk.start())) - pose.xyz(NamedAtomID("CA", prev.stop()));
	} else if ( has_next ) {
		const Loop& next = chunks_[i + 1];
		*bias = pose.xyz(NamedAtomID("CA", chunk.stop())) - pose.xyz(NamedAtomID("CA", next.start()));
	}

	// Normalize the vector and scale by the chainbreak bias term
	bias->normalize();
	*bias *= chainbreak_bias();
}

void RigidBodyMotionMover::update_chunks() {
	chunks_.clear();

	Size start = 1;
	for ( Size i = 1; i <= tree_.num_cutpoint(); ++i ) {
		Size stop = tree_.cutpoint(i);
		chunks_.add_loop(protocols::loops::Loop(start, stop));
		start = stop + 1;
	}

	chunks_.sequential_order();
}

double RigidBodyMotionMover::magnitude_rotation() const {
	return mag_rot_;
}

double RigidBodyMotionMover::magnitude_translation() const {
	return mag_trans_;
}

double RigidBodyMotionMover::chainbreak_bias() const {
	return cb_bias_;
}

const core::kinematics::FoldTree& RigidBodyMotionMover::fold_tree() const {
	return tree_;
}

void RigidBodyMotionMover::set_magnitude_rotation(double mag_rot) {
	debug_assert(mag_rot >= 0);
	mag_rot_ = mag_rot;
}

void RigidBodyMotionMover::set_magnitude_translation(double mag_trans) {
	debug_assert(mag_trans >= 0);
	mag_trans_ = mag_trans;
}

void RigidBodyMotionMover::set_chainbreak_bias(double cb_bias) {
	debug_assert(cb_bias >= 0);
	debug_assert(cb_bias <= 1);
	cb_bias_ = cb_bias;
}

void RigidBodyMotionMover::set_fold_tree(const core::kinematics::FoldTree& tree) {
	debug_assert(tree.check_fold_tree());
	tree_ = tree;
	update_chunks();
}

std::string RigidBodyMotionMover::get_name() const {
	return "RigidBodyMotionMover";
}

}  // namespace rigid
}  // namespace protocols
