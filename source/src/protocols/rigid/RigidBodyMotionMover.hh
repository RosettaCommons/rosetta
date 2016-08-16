// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rigid/RigidBodyMotionMover.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_RIGID_RIGIDBODYMOTIONMOVER_HH
#define INCLUDED_PROTOCOLS_RIGID_RIGIDBODYMOTIONMOVER_HH

// Unit header
#include <protocols/rigid/RigidBodyMotionMover.fwd.hh>

// C/C++ headers
#include <string>

// Project headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loops.hh>

// Package headers
#include <protocols/moves/Mover.hh>

//Auto Headers
#include <utility/vector1.hh>
namespace protocols {
namespace rigid {

class RigidBodyMotionMover : public moves::Mover {
public:
	explicit RigidBodyMotionMover(const core::kinematics::FoldTree& tree);

	/// @brief Randomly selects a chunk and perturbs its rigid body transformation
	void apply(core::pose::Pose& pose);

	/// @brief Returns the magnitude of translation
	double magnitude_translation() const;

	/// @brief Returns the magnitude of rotation
	double magnitude_rotation() const;

	/// @brief Returns the strength of the chainbreak bias
	double chainbreak_bias() const;

	/// @brief Returns the input fold tree
	const core::kinematics::FoldTree& fold_tree() const;

	/// @brief Updates the magnitude of translation
	void set_magnitude_translation(double mag_trans);

	/// @brief Updates the magnitude of rotation
	void set_magnitude_rotation(double mag_rot);

	/// @brief Updates the strength of the chainbreak bias
	void set_chainbreak_bias(double cb_bias);

	/// @brief Updates the input fold tree and regenerates chunks
	void set_fold_tree(const core::kinematics::FoldTree& tree);

	/// @brief Returns the name of this mover
	std::string get_name() const;

private:
	/// @brief Derives chunk definitions from input fold tree
	void update_chunks();

	/// @brief Having selected chunk i, compute the chainbreak-biased translation vector
	void compute_bias(unsigned i, const core::pose::Pose& pose, numeric::xyzVector<double>* cb_deltas) const;

	/// @brief Input fold tree
	core::kinematics::FoldTree tree_;

	/// @brief Chunks derived from input fold tree. Stored in increasing order of start position.
	protocols::loops::Loops chunks_;

	/// @brief Magnitude of rotation. Unless otherwise specified, equal to -rigid:rotation.
	double mag_rot_;

	/// @brief Magnitude of translation. Unless otherwise specified, equal to -rigid:translation.
	double mag_trans_;

	/// @brief Value on [0..1] that controls the strength of the chainbreak bias in translation moves
	double cb_bias_;
};

}  // namespace rigid
}  // namespace protocols

#endif  // PROTOCOLS_RIGID_RIGID_BODY_MOTION_MOVER_HH_
