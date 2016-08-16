// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/scmin/SidechainStateAssignment.hh
/// @brief  Declaration for the class that holds a state assignment for an entire system under continuous sidechain optimization
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_scmin_SidechainStateAssignment_HH
#define INCLUDED_core_pack_scmin_SidechainStateAssignment_HH

// Unit headers
#include <core/pack/scmin/SidechainStateAssignment.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/pack/scmin/AtomTreeCollection.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace core {
namespace pack {
namespace scmin {

/// @brief A simple class for tracking a network state and its energy where
/// each sidechain's state is described by a series of chi angles.
class SidechainStateAssignment : public utility::pointer::ReferenceCount {
public:
	SidechainStateAssignment( Size nmoltenres );
	SidechainStateAssignment( SidechainStateAssignment const & );
	SidechainStateAssignment &
	operator = ( SidechainStateAssignment const & );

	Size nmoltenres() const { return nmoltenres_; }
	scmin::ResidueAtomTreeCollectionMomento & state_momento( Size moltenresid );
	void assign_state( Size moltenresid, Size orig_rotid );
	void assign_energy( Real energy );

	scmin::ResidueAtomTreeCollectionMomento const & momento_for_moltenres( Size moltenresid ) const;
	Size orig_rotamer_id_for_moltenres( Size moltenres ) const;
	Real energy() const { return energy_; }
	bool any_unassigned() const;
	Size n_unassigned() const { return n_unassigned_; }

private:
	Size nmoltenres_;
	utility::vector1< scmin::ResidueAtomTreeCollectionMomento > state_assignments_;
	utility::vector1< Size > original_rotamer_id_;
	Real energy_;
	Size n_unassigned_;
};


} // namespace scmin
} // namespace pack
} // namespace core

#endif
