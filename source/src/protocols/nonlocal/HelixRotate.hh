// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/HelixRotate.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_protocols_nonlocal_HelixRotate_HH
#define INCLUDED_protocols_nonlocal_HelixRotate_HH

// Unit header
#include <protocols/nonlocal/HelixRotate.fwd.hh>

// C/C++ headers
#include <string>

// Utility headers

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loop.hh>

// Package headers
#include <protocols/moves/Mover.hh>

//Auto Headers
#include <protocols/loops/Loops.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace nonlocal {

/// @class A simple protocol for rotating a contiguous stretch of residues along the
/// axis defined by its first and last residues. Makes no assumptions about the
/// kinematics of the system.
///
/// TODO(cmiles) Improved handling of kinked helices
class HelixRotate : public moves::Mover {
public:
	HelixRotate();
	HelixRotate(const protocols::loops::Loop& helix, double degrees);

	/// @brief Rotates the helix by the specified number of degrees
	void apply(core::pose::Pose& pose);

	/// @brief Returns the protocol's name
	std::string get_name() const;

	/// @brief Creates a new instance with the copy constructor
	moves::MoverOP clone() const;

	/// @brief Creates a new instance with the default constructor
	moves::MoverOP fresh_instance() const;

	/// @brief Returns the helix to be modified
	const protocols::loops::Loop& get_helix() const;

	/// @brief Returns the number of degrees to rotate the helix
	double get_degrees() const;

	/// @brief Updates the helix to be modified
	void set_helix(const protocols::loops::Loop& helix);

	/// @brief Updates the number of degrees to rotate the helix
	void set_degrees(double degrees);

private:
	/// @brief Shared initialization routine
	void initialize(const protocols::loops::Loop& helix, double degrees);

	/// @brief Returns true if this instance is valid and fully configured
	bool is_valid() const;

	/// @brief Partitions the structure such that the helix to be rotated belongs
	/// to its own chunk, which will be subsequently connected to the star fold
	/// tree by its own jump
	void decompose_structure(unsigned num_residues, protocols::loops::Loops* chunks) const;

	/// @brief Searches chunks for the member representing the helix, returning its index
	unsigned jump_containing_helix(const protocols::loops::Loops& chunks) const;

	/// @brief Computes rotational parameters-- axis and point
	void get_rotation_parameters(
		const core::pose::Pose& pose,
		numeric::xyzVector<double>* axis,
		numeric::xyzVector<double>* point) const;

	/// @brief Stretch of contiguous residues representing the helix to be rotated
	protocols::loops::Loop helix_;

	/// @brief Number of degrees to rotate the helix
	double degrees_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif
