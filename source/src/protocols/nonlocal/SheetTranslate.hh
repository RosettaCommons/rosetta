// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/SheetTranslate.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_protocols_nonlocal_SheetTranslate_HH
#define INCLUDED_protocols_nonlocal_SheetTranslate_HH

// Unit header
#include <protocols/nonlocal/SheetTranslate.fwd.hh>

// C/C++ headers
#include <string>

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

/// @class A simple protocol for translating a contiguous stretch of residues along
/// the axis defined by its first and last residues. Makes no assumptions about the
/// kinematics of the system.
///
/// TODO(cmiles) Consider having the user define the axis of translation.
/// TODO(cmiles) Improved handling of curved sheets
class SheetTranslate : public moves::Mover {
public:
	SheetTranslate();
	SheetTranslate(const protocols::loops::Loop& sheet, double distance_ang);

	/// @brief Translates the sheet by the specified distance (in Angstroms)
	void apply(core::pose::Pose& pose);

	/// @brief Returns the protocol's name
	std::string get_name() const;

	/// @brief Creates a new instance with the copy constructor
	moves::MoverOP clone() const;

	/// @brief Creates a new instance with the default constructor
	moves::MoverOP fresh_instance() const;

	/// @brief Returns the sheet to be modified
	const protocols::loops::Loop& get_sheet() const;

	/// @brief Returns the distance (in Angstroms) to translate the sheet
	double get_distance() const;

	/// @brief Updates the sheet to be modified
	void set_sheet(const protocols::loops::Loop& sheet);

	/// @brief Updates the distance (in Angstroms) to translate the sheet
	void set_distance(double distance_ang);

private:
	/// @brief Shared initialization routine
	void initialize(const protocols::loops::Loop& sheet, double distance);

	/// @brief Returns true if this instance is valid and fully configured
	bool is_valid() const;

	/// @brief Partitions the structure such that the sheet to be translated
	/// belongs to its own chunk, which will be subsequently connected to the
	/// star fold tree by its own jump
	void decompose_structure(unsigned num_residues, protocols::loops::Loops* chunks) const;

	/// @brief Searches chunks for the member representing the sheet, returning its index
	unsigned jump_containing_sheet(const protocols::loops::Loops& chunks) const;

	/// @brief Stretch of contiguous residues representing the sheet to be rotated
	protocols::loops::Loop sheet_;

	/// @brief Distance to translate the sheet in Angstroms
	double distance_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif
