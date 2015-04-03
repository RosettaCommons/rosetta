// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamer_set/RotamerSubset.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/pack/rotamer_set/RotamerSubset.hh>

// Project Headers
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/methods/Methods.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <iostream>

#include <core/conformation/AbstractRotamerTrie.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace rotamer_set {

static thread_local basic::Tracer tt( "core.pack.rotamer_set.RotamerSubset", basic::t_info );

RotamerSubset::RotamerSubset(
	RotamerSet & src_rotset,
	utility::vector1< Size > const & rotamer_subset
) :
	n_residue_types_( 0 ),
	cached_tries_( scoring::methods::n_energy_methods, 0 ),
	id_for_current_rotamer_( 0 ),
	rotamer_offsets_require_update_( false )
{
	rotamers_.reserve( rotamer_subset.size() );
	Size const src_curr_rot_index = src_rotset.id_for_current_rotamer();
	for ( Size ii = 1; ii <= rotamer_subset.size(); ++ii ) {
		steal_rotamer( src_rotset.nonconst_rotamer( rotamer_subset[ ii ] ) );
		if ( rotamer_subset[ ii ] == src_curr_rot_index ) {
			id_for_current_rotamer_ = rotamers_.size();
		}
	}

}

RotamerSubset::~RotamerSubset() {}

void
RotamerSubset::add_rotamer(
	conformation::Residue const & rotamer
)
{
	prepare_for_new_residue_type( rotamer.type() );
	push_back_rotamer( rotamer.clone() );
}


Size
RotamerSubset::get_n_residue_types() const
{
	update_rotamer_offsets();
	return n_residue_types_;
}

Size
RotamerSubset::get_n_residue_groups() const
{
	update_rotamer_offsets();
	return n_residue_groups_;
}


Size
RotamerSubset::get_residue_type_begin( Size which_restype ) const
{
	update_rotamer_offsets();
debug_assert( which_restype <= n_residue_types_ );
	return residue_type_rotamers_begin_[ which_restype ];
}

Size
RotamerSubset::get_residue_group_begin( Size which_resgroup ) const
{
	update_rotamer_offsets();
debug_assert( which_resgroup <= n_residue_groups_ );
	return residue_group_rotamers_begin_[ which_resgroup ];
}


Size
RotamerSubset::get_n_rotamers_for_residue_type( Size which_restype ) const
{
	update_rotamer_offsets();

debug_assert( which_restype <= n_residue_types_ );
	return n_rotamers_for_restype_[ which_restype ];
}

Size
RotamerSubset::get_n_rotamers_for_residue_group( Size which_resgroup ) const
{
	update_rotamer_offsets();

debug_assert( which_resgroup <= n_residue_groups_ );
	return n_rotamers_for_resgroup_[ which_resgroup ];
}

/// @brief given a rotamer id, return an int which represents a type for this rotamer.
/// INCOMPLETELY IMPLEMENTED. ANDREW: FIX THIS.
Size
RotamerSubset::get_residue_type_index_for_rotamer( Size /* which_rotamer */ ) const
{
	//return residue_type_for_rotamers_[ which_rotamer ];
	return 1;
}

/// @brief given a rotamer id, return an int which represents a type for this rotamer.
/// INCOMPLETELY IMPLEMENTED. ANDREW: FIX THIS.
Size
RotamerSubset::get_residue_group_index_for_rotamer( Size /* which_rotamer */ ) const
{
	//return residue_type_for_rotamers_[ which_rotamer ];
	return 1;
}


Size
RotamerSubset::num_rotamers() const
{
	return rotamers_.size();
}

Size
RotamerSubset::id_for_current_rotamer() const
{
	return id_for_current_rotamer_;
}

conformation::ResidueCOP
RotamerSubset::rotamer( Size rot_id ) const
{
	return rotamers_[ rot_id ];
}

conformation::Residue const &
RotamerSubset::rotamer_ref( Size rot_id ) const
{
	return *rotamers_[ rot_id ];
}


/// @details In handing out non-const data, the guarantee of rotamer-type contiguity
/// within the rotamers_ array, and the correspondence of the rotamer offset
/// data is lost.  Future access to rotamer offset data first requires an update
/// of the rotamer offset arrays.
conformation::ResidueOP
RotamerSubset::nonconst_rotamer( Size rot_id )
{
	rotamer_offsets_require_update_ = true;

	return rotamers_[ rot_id ];
}

void
RotamerSubset::store_trie(
	Size method_enum_id,
	conformation::AbstractRotamerTrieOP trie
)
{
	cached_tries_[ method_enum_id ] = trie;
}


conformation::AbstractRotamerTrieCOP
RotamerSubset::get_trie( Size method_enum_id ) const
{
	return cached_tries_[ method_enum_id ];
}

/// @details  O(n) operation; if you have a lot of rotamers you want to remove, use
/// drop_rotamers() instead.
void
RotamerSubset::drop_rotamer( Size rot_id )
{
debug_assert( rot_id <= rotamers_.size() );
	utility::vector1< conformation::ResidueOP > copy_rotamers( rotamers_.size() - 1, 0 );
	Size count_copy( 1 );
	for ( Size ii = 1; ii <= rotamers_.size(); ++ii ) {
		if ( ii != rot_id ) {
			copy_rotamers[ count_copy ] = rotamers_[ ii ];
			if ( ii == id_for_current_rotamer_ ) {
				id_for_current_rotamer_ = count_copy;
			}
			++count_copy;
		} else {
			if ( ii == id_for_current_rotamer_ ) {
				id_for_current_rotamer_ = 0;
			}
		}
	}
	copy_rotamers.swap( rotamers_ );
	rotamer_offsets_require_update_ = true;
	update_rotamer_offsets();

}

/// @brief rotamers_to_delete must be of size nrotmaers -- each position
/// in the array that's "true" is removed from the set of rotamers
void
RotamerSubset::drop_rotamers( utility::vector1< bool > const & rotamers_to_delete )
{
debug_assert( rotamers_to_delete.size() == rotamers_.size() );

	Size n_dropped = 0;
	for ( Size ii = 1; ii <= rotamers_.size(); ++ii ) {
		if ( rotamers_to_delete[ ii ] ) {
			/// if all rotamers end up dropped, then preserve the input rotamer.
			if ( ii == id_for_current_rotamer_ ) {
				current_rotamer_copy_ = rotamers_[ ii ];
			}
			rotamers_[ ii ] = 0;
			++n_dropped;
		}
	}
	if ( n_dropped == 0 ) return;

	if ( n_dropped == rotamers_.size() ) {
		if ( id_for_current_rotamer_ == 0 ) {
			utility_exit_with_message( "ERROR:: RotamerSubset::drop_rotamers attempted to remove all rotamers without available input_rotamer." );
		}
		// keep the input rotamer.
		rotamers_.resize( 1 );
		rotamers_[ 1 ] = current_rotamer_copy_;
		id_for_current_rotamer_ = 1;
		current_rotamer_copy_.reset();
	} else {
		utility::vector1< conformation::ResidueOP > new_rotamers( rotamers_to_delete.size() - n_dropped, 0 );
		Size count_new = 1;
		for ( Size ii = 1; ii <= rotamers_.size(); ++ii ) {
			if ( rotamers_[ ii ] != 0 ) {
				new_rotamers[ count_new ] = rotamers_[ ii ];
				if ( ii == id_for_current_rotamer_ ) {
					id_for_current_rotamer_ = count_new;
				}
				++count_new;
			}
		}
		new_rotamers.swap( rotamers_ );
	}
	rotamer_offsets_require_update_ = true;
	update_rotamer_offsets();
}

/// @brief deletes the rotamers in the list with the given indices.
/// The indices of these rotamers is presumed to be those before any delete operation.
/// e.g. if there are four rotamers, and rotamer_indices_to_delete includes 1 & 3,
/// then the rotamers that will remain are the rotamers originally indexed as 2 and 4,
/// even though their new indices will be 1 & 2.
void
RotamerSubset::drop_rotamers_by_index(
	utility::vector1< Size > const & rotamer_indices_to_delete
)
{
	utility::vector1< bool > rotamers_to_delete( rotamers_.size(), false );
	for ( Size ii = 1; ii <= rotamer_indices_to_delete.size(); ++ii ) {
		rotamers_to_delete[ rotamer_indices_to_delete[ ii ] ] = true;
	}
	drop_rotamers( rotamers_to_delete );
}

void
RotamerSubset::prepare_for_new_residue_type( core::chemical::ResidueType const & restype )
{
	if ( n_residue_types_ == 0 ) {
		new_residue_type();
		new_residue_group();
		return;
	}

	if ( different_restype( rotamers_[ num_rotamers() ]->type(), restype )) {
		new_residue_type();
	}
	if (  different_resgroup( rotamers_[ num_rotamers() ]->type(), restype )) {
		new_residue_group();
	}
}

bool
RotamerSubset::different_restype( core::chemical::ResidueType const & rt1, core::chemical::ResidueType const & rt2 ) const
{
	return & rt1 != & rt2;
}

/// @details The logic to determine if two residue types should be classified as part of the same group.
/// The thinking is as follows.  Two residue types are in the same group if they have the same residue type.
/// They're in the same group if their residue types differ, but they have the same name3 (HIS vs HIS_D have
/// the same name3) and they have the same neighbor radius (SER and PhosphoSER should have different groups).
/// The goal is to organize residue types together which will be packed together (as happens in multistate design
/// with HIS and HISD) and that have the same reach (as is needed for the AANeighborSparseMatrix).
bool
RotamerSubset::different_resgroup( core::chemical::ResidueType const & rt1, core::chemical::ResidueType const & rt2 ) const
{
	return & rt1 != & rt2 && ( rt1.name3() != rt2.name3() || rt1.nbr_radius() != rt2.nbr_radius() );
}

void
RotamerSubset::new_residue_type()
{
	++n_residue_types_;
	residue_type_rotamers_begin_.push_back( num_rotamers() + 1);
	n_rotamers_for_restype_.push_back( 0 );
}

void
RotamerSubset::new_residue_group()
{
	++n_residue_groups_;
	residue_group_rotamers_begin_.push_back( num_rotamers() + 1 );
	n_rotamers_for_resgroup_.push_back( 0 );
}

void
RotamerSubset::push_back_rotamer( conformation::ResidueOP rotamer )
{
	rotamers_.push_back( rotamer );
	residue_type_for_rotamers_.push_back( n_residue_types_ );
	++n_rotamers_for_restype_[ n_residue_types_ ];
}

void
RotamerSubset::update_rotamer_offsets() const
{
	if ( ! rotamer_offsets_require_update_ ) return;

	if ( rotamers_.size() == 0 ) {
		n_residue_types_ = 0;
		n_residue_groups_ = 0;
		residue_type_for_rotamers_.resize( 0 );
		residue_group_for_rotamers_.resize( 0 );
		residue_type_rotamers_begin_.resize( 0 );
		residue_group_rotamers_begin_.resize( 0 );
		n_rotamers_for_restype_.resize( 0 );
		n_rotamers_for_resgroup_.resize( 0 );
		return;
	}

	/// From here forward, rotamers_.size() >= 1
	residue_type_for_rotamers_.resize( rotamers_.size() );
	residue_group_for_rotamers_.resize( rotamers_.size() );
	n_residue_types_ = 1;
	n_residue_groups_ = 1;
	residue_type_for_rotamers_[ 1 ] = n_residue_types_;
	residue_group_for_rotamers_[ 1 ] = n_residue_groups_;
	for ( Size ii = 2; ii <= rotamers_.size(); ++ii ) {
		// compare addresses of the two types
		// treat them as different amino acids only if they have different name3's
		// or if they have different radii
		//if ( & (rotamers_[ ii ]->type()) != & (rotamers_[ ii ]->type()) ) {
		if ( different_restype( rotamers_[ ii ]->type(), rotamers_[ ii-1 ]->type() ) ) {
			++n_residue_types_;
		}
		residue_type_for_rotamers_[ ii ] = n_residue_types_;

		if ( different_resgroup( rotamers_[ ii ]->type(), rotamers_[ ii-1 ]->type() ) ) {
			++n_residue_groups_;
		}
		residue_group_for_rotamers_[ ii ] = n_residue_groups_;
	}

	residue_type_rotamers_begin_.resize( n_residue_types_ );
	n_rotamers_for_restype_.resize( n_residue_types_ );
	std::fill( residue_type_rotamers_begin_.begin(), residue_type_rotamers_begin_.end(), 0 );
	std::fill( n_rotamers_for_restype_.begin(), n_rotamers_for_restype_.end(), 0 );

	residue_group_rotamers_begin_.resize( n_residue_groups_ );
	n_rotamers_for_resgroup_.resize( n_residue_groups_ );
	std::fill( residue_group_rotamers_begin_.begin(), residue_group_rotamers_begin_.end(), 0 );
	std::fill( n_rotamers_for_resgroup_.begin(), n_rotamers_for_resgroup_.end(), 0 );

	Size count_seen_residue_types( 1 );
	Size count_seen_residue_groups( 1 );
	n_rotamers_for_restype_[ count_seen_residue_types ] = 1;
	n_rotamers_for_resgroup_[ count_seen_residue_groups ] = 1;
	residue_type_rotamers_begin_[ count_seen_residue_types ] = 1;
	residue_group_rotamers_begin_[ count_seen_residue_groups ] = 1;

	for ( Size ii = 2; ii <= rotamers_.size(); ++ii ) {
		if ( residue_type_for_rotamers_[ ii ] != residue_type_for_rotamers_[ ii-1 ] ) {
			++count_seen_residue_types;
			residue_type_rotamers_begin_[ count_seen_residue_types ] = ii;
		}
		++n_rotamers_for_restype_[ count_seen_residue_types ];
		if ( residue_group_for_rotamers_[ ii ] != residue_group_for_rotamers_[ ii-1 ] ) {
			++count_seen_residue_groups;
			residue_group_rotamers_begin_[ count_seen_residue_groups ] = ii;
		}
		++n_rotamers_for_resgroup_[ count_seen_residue_groups ];
	}
	//std::cout << "nrestypes " << n_residue_types_ << std::endl;
	rotamer_offsets_require_update_ = false;
}

void
RotamerSubset::build_rotamers(
	pose::Pose const & ,
	scoring::ScoreFunction const &,
	task::PackerTask const &,
	graph::GraphCOP,
	bool
)
{}

void
RotamerSubset::build_dependent_rotamers(
	RotamerSets const &,
	pose::Pose const &,
	scoring::ScoreFunction const &,
	task::PackerTask const &,
	graph::GraphCOP
)
{}

void
RotamerSubset::compute_one_body_energies(
	pose::Pose const &,
	scoring::ScoreFunction const &,
	task::PackerTask const &,
	graph::GraphCOP,
	utility::vector1< core::PackerEnergy > &
) const
{}

void
RotamerSubset::compute_one_and_two_body_energies(
	pose::Pose const &,
	scoring::ScoreFunction const &,
	task::PackerTask const &,
	graph::GraphCOP,
	utility::vector1< core::PackerEnergy > &,
	utility::vector1< utility::vector1< core::PackerEnergy > > &,
	utility::vector1< core::Size > &
) const
{}

void
RotamerSubset::compute_one_body_energy_maps(
	pose::Pose const &,
	scoring::ScoreFunction const &,
	task::PackerTask const &,
	graph::GraphCOP,
	utility::vector1< scoring::EnergyMap > &
) const
{}


/// @brief declare that a new block of residue types has begun, and that new residues
/// are about to be pushed back.
//void
//RotamerSubset::declare_new_residue_type()
//{
//	++n_residue_types_;
//	residue_type_rotamers_begin_.push_back( num_rotamers() + 1);
//	n_rotamers_for_restype_.push_back( 0 );
//}

/// @brief appends a rotamer to the list of rotamers, and increments the count
/// for the number of rotamers for the current value of n_residue_types.
//void
//RotamerSubset::push_back_rotamer( conformation::ResidueOP rotamer )
//{
//	rotamers_.push_back( rotamer );
//	residue_type_for_rotamers_.push_back( n_residue_types_ );
//	++n_rotamers_for_restype_[ n_residue_types_ ];
//}


void RotamerSubset::steal_rotamer( conformation::ResidueOP rotamer ) {
	rotamer_offsets_require_update_ = true;
	rotamers_.push_back( rotamer );
}


} // rotamer_set
} // pack
} // core
