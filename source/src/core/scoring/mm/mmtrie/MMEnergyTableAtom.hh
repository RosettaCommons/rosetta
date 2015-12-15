// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/mmtrie/MMEnergyTableAtom.hh
/// @brief  Header for the MMEnergyTableAtom. Heavily coppied from the EtableAtom.hh
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

#ifndef INCLUDED_core_scoring_mm_mmtrie_MMEnergyTableAtom_hh
#define INCLUDED_core_scoring_mm_mmtrie_MMEnergyTableAtom_hh

// Unit Headers
#include <core/scoring/mm/mmtrie/MMEnergyTableAtom.fwd.hh>

// Project Headers
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/types.hh>

// STL Headers
#include <iostream>

// Numceric Headers
#include <numeric/xyzVector.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace mm {
namespace mmtrie {

class MMEnergyTableAtom : public conformation::Atom
{
public:
	typedef conformation::Atom parent;

public:
	MMEnergyTableAtom();

	MMEnergyTableAtom( conformation::Residue const & res, Size mm_atom_index );

	/// @brief non-virtual destructor to keep MMEnergyTableAtom small and lightweight
	/// as a virtual destructor would add a vtable to the class
	/// But I fear leaks... do I know how xyzVector dealloates its data?
	virtual ~MMEnergyTableAtom();

	/// @brief deprecated!
	inline
	Size mm_atom_type() const { return mm_type(); }
	/// @brief deprecated!
	inline
	void mm_atom_type( Size setting ) { mm_type( setting ); }

	/// @brief property required by RotamerTrie class
	inline
	bool is_hydrogen() const { return is_hydrogen_; }

	/// @brief setter method for data required by RotamerTrie class
	inline
	void is_hydrogen( bool setting ) { is_hydrogen_ = setting; }

	/// @brief send a description of the atom to standard out
	void print() const;

	/// @brief send a description of the atom to an output stream
	void print( std::ostream & os ) const;

	/// @brief compairison operator for sorting
	inline
	bool operator < ( MMEnergyTableAtom const & other ) const
	{
		if ( mm_type() == other.mm_type() ) {
			// Quoting Jack Snoeyink: "Epsilons in code are evil."  But whatcha gonna do?
			// In this case, though, you could get points a, b, and c such that a == b, b == c, but a < c.
			// In rare cases this would cause the std::sort() insertion sort step
			// to run off the end of the array/vector and cause a segfault in the trie.
			for ( int ii = 0; ii< 3; ++ii ) {
				//Distance diff = std::abs( xyz()[ ii ] - other.xyz()[ ii ] );
				//if ( diff > EPSILON ) {
				if ( float(xyz()[ ii ]) != float(other.xyz()[ ii ]) ) {
					return xyz()[ ii ] < other.xyz()[ ii ];
				}
			}
		} else {
			return mm_type() < other.mm_type();
		}
		return false;
	}

	/// @brief equality operator for shared-prefix detection
	inline
	bool operator == ( MMEnergyTableAtom const & other ) const
	{
		if ( mm_type() == other.mm_type() ) {
			for ( int ii = 0; ii< 3; ++ii ) {
				// Epsilons bad -- see above.
				//Distance diff = std::abs( xyz()[ ii ] - other.xyz()[ ii ] );
				//if ( diff > EPSILON ) {
				if ( float(xyz()[ ii ]) != float(other.xyz()[ ii ]) ) {
					return false;
				}
			}
		} else {
			return false;
		}
		return true;
	}


private:

	bool is_hydrogen_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

std::ostream & operator << ( std::ostream & os, MMEnergyTableAtom const & atom );

} // namespace mmtrie
} // namespace mm
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_mm_mmtrie_MMEnergyTableAtom )
#endif // SERIALIZATION


#endif
