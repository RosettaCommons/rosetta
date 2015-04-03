// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/ResidueCoordinateChangeList.fwd.hh
/// @brief  AtomTree/Conformation communication vector class forward declaration
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_kinematics_ResidueCoordinateChangeList_hh
#define INCLUDED_core_kinematics_ResidueCoordinateChangeList_hh

// Unit headers
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>

// Package headers
#include <core/id/AtomID.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>


namespace core {
namespace kinematics {

/// @brief The AtomTree is responsible for informing the conformation
/// of which residues have had either internal (DOF) or external
/// (xyz) coordinate changes so that the Conformation may shuttle
/// O(k) -- output sensitive -- data from the AtomTree to the
/// Residue objects it manages.
class ResidueCoordinateChangeList : public utility::pointer::ReferenceCount
{
public:
	ResidueCoordinateChangeList();

	virtual ~ResidueCoordinateChangeList();

	ResidueCoordinateChangeList const &
	operator = ( ResidueCoordinateChangeList const & rhs );

	/// @brief Set the number of residues in the conformation being tracked.
	void
	total_residue( Size total_residue );

	/// @brief Return the number of
	Size
	total_residue() const {
		return total_residue_;
	}

	/// @brief Is the list of the changed residues empty?
	bool
	empty() const;

	/// @brief Reset the list of those residues that have moved.  O(k).
	void
	clear();

	/// @brief Mark a residue as having changed by passing in an AtomId for one atom
	/// in that residue
	void
	mark_residue_moved( id::AtomID atid );

	/// @brief Mark a residue as having changed by passing in the index of that residue.
	void
	mark_residue_moved( Size resid );

	/// @brief returns an iterator to the beginning of the (unordered) list of
	/// residues that have moved.
	ResidueListIterator
	residues_moved_begin() const;

	/// @brief returns an iterator to the end of the (unordered) list of
	/// residues that have moved.
	ResidueListIterator
	residues_moved_end() const;

private:
	/// @brief O(N) -- used in assert statements in debug mode
	bool
	no_residues_with_nonzero_change_id() const;

private:
	ResidueIndexList changed_residues_;

	/// @brief The reverse mapping from residue index to its position in the
	/// list of residues that have changed.  Size total_residue_
	utility::vector1< Size > residue_change_id_;

	Size total_residue_;

};


} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_ResidueCoordinateChangeList_HH
