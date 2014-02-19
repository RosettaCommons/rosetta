// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/scmin/SidechainStateAssignment.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/pack/scmin/SidechainStateAssignment.hh>

// Package Headers
#include <core/pack/scmin/AtomTreeCollection.hh>

// AUTO-REMOVED #include <core/conformation/Residue.hh> // REQUIRED FOR WINDOWS

#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace scmin {

SidechainStateAssignment::SidechainStateAssignment( Size nmoltenres ) :
	nmoltenres_( nmoltenres ),
	state_assignments_( nmoltenres ),
	original_rotamer_id_( nmoltenres, 0 ),
	energy_( 0.0 ),
	n_unassigned_( nmoltenres )
{}

SidechainStateAssignment::SidechainStateAssignment( SidechainStateAssignment const & src ) :
		utility::pointer::ReferenceCount(src),
		nmoltenres_( src.nmoltenres_ ),
		state_assignments_( src.state_assignments_ ),
		original_rotamer_id_( src.original_rotamer_id_ ),
		energy_( src.energy_ ),
		n_unassigned_( src.n_unassigned_ )
{}

SidechainStateAssignment const &
SidechainStateAssignment::operator = ( SidechainStateAssignment const & rhs )
{
	if ( this != &rhs ) {
		nmoltenres_ =  rhs.nmoltenres_;
		state_assignments_ = rhs.state_assignments_;
		original_rotamer_id_ = rhs.original_rotamer_id_;
		energy_ = rhs.energy_;
		n_unassigned_ = rhs.n_unassigned_;
	}
	return *this;
}

scmin::ResidueAtomTreeCollectionMomento &
SidechainStateAssignment::state_momento( Size moltenresid )
{
	return state_assignments_[ moltenresid ];
}

void
SidechainStateAssignment::assign_state(
	Size moltenresid,
	Size orig_rotid
)
{
	assert( orig_rotid != 0 );
	if ( original_rotamer_id_[ moltenresid ] == 0 ) {
		assert( n_unassigned_ > 0 );
		--n_unassigned_;
	}
	original_rotamer_id_[ moltenresid ] = orig_rotid;
}

void SidechainStateAssignment::assign_energy( Real energy ) { energy_ = energy; }

scmin::ResidueAtomTreeCollectionMomento const &
SidechainStateAssignment::momento_for_moltenres( Size moltenresid ) const
{
	return state_assignments_[ moltenresid ];
}


Size
SidechainStateAssignment::orig_rotamer_id_for_moltenres( Size moltenres ) const
{
	return original_rotamer_id_[ moltenres ];
}

bool
SidechainStateAssignment::any_unassigned() const
{
	return n_unassigned_ != 0;
}

} // namespace scmin
} // namespace pack
} // namespace core

