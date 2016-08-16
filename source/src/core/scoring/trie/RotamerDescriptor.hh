// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/trie/RotamerDescriptor.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_trie_RotamerDescriptor_hh
#define INCLUDED_core_scoring_trie_RotamerDescriptor_hh

#include <core/types.hh>


#include <utility/vector1.fwd.hh>

namespace core {
namespace scoring {
namespace trie {

template < class AT, class CPDAT >
class RotamerDescriptorAtom
{
public:
	RotamerDescriptorAtom() {}
	RotamerDescriptorAtom( AT atom, CPDAT cp_data ) : atom_( atom ), cp_data_( cp_data ) {}

	AT const & atom() const { return atom_; }
	CPDAT const & cp_data() const { return cp_data_; }

	void atom( AT setting ) { atom_ = setting; }
	void cp_data( CPDAT setting ) { cp_data_ = setting; }

	bool operator < ( RotamerDescriptorAtom< AT, CPDAT > const & other ) const
	{
		if ( atom_ < other.atom_ ) return true;
		else if ( atom_ == other.atom_ && cp_data_ < other.cp_data_ ) return true;
		return false;
	}

	bool operator == ( RotamerDescriptorAtom< AT, CPDAT > const & other ) const
	{
		return (atom_ == other.atom_) && (cp_data_ == other.cp_data_ );
	}

private:
	AT atom_;
	CPDAT cp_data_;
};


template < class AT, class CPDAT >
class RotamerDescriptor
{
public:
	RotamerDescriptor() : rotamer_id_( 0 ), natoms_( 0 ) {}

	// Setters
	void natoms( Size setting ) { natoms_ = setting; atoms_.resize( setting ); }
	void atom( Size index, RotamerDescriptorAtom< AT, CPDAT > const & newatom ) { atoms_[ index ] = newatom; }
	void rotamer_id( Size setting ) { rotamer_id_ = setting; }

	// Getters
	Size natoms() const { return natoms_; }
	Size rotamer_id() const { return rotamer_id_; }
	RotamerDescriptorAtom< AT, CPDAT > const & atom( Size index ) { return atoms_[ index ]; }

	// Comparison operator for sorting.
	bool operator < ( RotamerDescriptor< AT, CPDAT > const & other ) const
	{
		if ( natoms_ < other.natoms_ ) { return true; }
		else if ( natoms_ > other.natoms_ ) { return false; }
		else {
			for ( Size ii = 1; ii <= natoms_; ++ii ) {
				if ( atoms_[ ii ] < other.atoms_[ ii ] ) {
					//std::cout << "atom " << ii << "lt" << std::endl;
					//atoms_[ ii ].atom().print();
					//other.atoms_[ ii ].atom().print();
					return true;
				} else if ( other.atoms_[ ii ] < atoms_[ ii ] ) {
					//std::cout << "atom " << ii << "gt" << std::endl;
					//atoms_[ ii ].atom().print();
					//other.atoms_[ ii ].atom().print();
					return false;
				}
			}
		}
		return false;
	}


	Size
	count_atoms_in_common( RotamerDescriptor< AT, CPDAT > const & other ) const
	{
		Size const n_to_compare = natoms_ < other.natoms_ ? natoms_ : other.natoms_;
		for ( Size ii = 1; ii <= n_to_compare; ++ii ) {
			if ( ! (atoms_[ ii ] == other.atoms_[ ii ]) ) {
				return ii - 1;
			}
		}
		return n_to_compare;
	}

private:
	Size rotamer_id_;
	typename utility::vector1< RotamerDescriptorAtom< AT, CPDAT > > atoms_;
	Size natoms_;

};


} // namespace trie
} // namespace scoring
} // namespace core

#endif

