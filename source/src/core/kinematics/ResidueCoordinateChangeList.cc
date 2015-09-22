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

// Unit headers
#include <core/kinematics/ResidueCoordinateChangeList.hh>

// Package headers
#include <core/id/AtomID.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer tr( "core.kinematics", basic::t_info );

namespace core {
namespace kinematics {

ResidueCoordinateChangeList::ResidueCoordinateChangeList()
:
	ReferenceCount(),
	total_residue_( 0 )
{
}

ResidueCoordinateChangeList::~ResidueCoordinateChangeList()
{}

ResidueCoordinateChangeList const &
ResidueCoordinateChangeList::operator = ( ResidueCoordinateChangeList const & rhs )
{
	if ( this == &rhs ) return *this;

	clear(); // O(k) -- wipe change information from up until now.

	total_residue_ = rhs.total_residue_;
	changed_residues_.reserve( total_residue_ );
	changed_residues_ = rhs.changed_residues_;
	residue_change_id_.resize( total_residue_, 0 );
	if ( changed_residues_.size() != 0 ) {
		// O(k) -- copy over the parts that need changing
		for ( Size ii = 1; ii <= changed_residues_.size(); ++ii ) {
			residue_change_id_[ changed_residues_[ ii ] ] = rhs.residue_change_id_[ changed_residues_ [ ii ] ];
		}
	}
	return *this;
}


/// @details adding and deleting residues repeatedly should be O(1)
/// so use the stl resize-with-argument command when resizing residue_change_id_.
/// That way, if residue_change_id_ has reserved space of at least total_residue,
/// then the resize will initialize exactly the number of new residues that have
/// been added.
/// In debug mode, this operation is O( total_residue ).
void
ResidueCoordinateChangeList::total_residue( Size total_residue )
{
	debug_assert( empty() );
	debug_assert( no_residues_with_nonzero_change_id() );

	total_residue_ = total_residue;
	changed_residues_.reserve( total_residue_ );
	residue_change_id_.resize( total_residue_, 0 );
}


/// @brief Is the list of the changed residues empty? bool
bool
ResidueCoordinateChangeList::empty() const
{
	return changed_residues_.size() == 0;
}

bool
ResidueCoordinateChangeList::no_residues_with_nonzero_change_id() const
{
	for ( Size ii = 1; ii <= total_residue_; ++ii ) {
		if ( residue_change_id_[ ii ] != 0 ) {
			return false;
		}
	}
	return true;
}

/// @details -- reset the residue_change_id_ for those residues which had previously changed
void
ResidueCoordinateChangeList::clear()
{
	for ( Size ii = 1; ii <= changed_residues_.size(); ++ii ) {
		debug_assert( residue_change_id_[ changed_residues_[ ii ] ] != 0 );
		residue_change_id_[ changed_residues_[ ii ] ] = 0;
	}
	changed_residues_.clear();
}

/// @brief Mark a residue as having changed by passing in an AtomId for one atom
/// in that residue
void
ResidueCoordinateChangeList::mark_residue_moved( id::AtomID atid )
{
	mark_residue_moved( atid.rsd() );
}

/// @brief Mark a residue as having changed by passing in the index of that residue.
void
ResidueCoordinateChangeList::mark_residue_moved( Size resid )
{
	if ( !(resid>0 && resid <= total_residue_ ) ) {
		tr.Error << "ASSERTION FAILED IN ResidueCoordinateChangeList: resid: " << resid << " total_residue " << total_residue_ << std::endl;
	}
	debug_assert( resid > 0 && resid <= total_residue_ );
	if ( residue_change_id_[ resid ] == 0 ) {
		changed_residues_.push_back( resid );
		residue_change_id_[ resid ] = changed_residues_.size();
	} /// else, the residue has already been marked as having changed

}

/// @brief returns an iterator to the beginning of the (unordered) list of
/// residues that have moved.
ResidueListIterator
ResidueCoordinateChangeList::residues_moved_begin() const
{
	return changed_residues_.begin();
}

/// @brief returns an iterator to the end of the (unordered) list of
/// residues that have moved.
ResidueListIterator
ResidueCoordinateChangeList::residues_moved_end() const
{
	return changed_residues_.end();
}


} // namespace kinematics
} // namespace core

