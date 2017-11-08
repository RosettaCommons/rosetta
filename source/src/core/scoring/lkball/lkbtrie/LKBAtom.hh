// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/etable/etrie/EtableAtom.hh
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_lkball_lkbtrie_LKBAtom_hh
#define INCLUDED_core_scoring_lkball_lkbtrie_LKBAtom_hh

// Unit Headers
#include <core/scoring/lkball/lkbtrie/LKBAtom.fwd.hh>

// Package headres
#include <core/scoring/lkball/LK_BallInfo.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Atom.hh>

// STL Headers
#include <iostream>

// Numceric Headers
#include <numeric/xyzVector.hh>

namespace core {
namespace scoring {
namespace lkball {
namespace lkbtrie {

class LKBAtom
{
public:
	LKBAtom();

	~LKBAtom();

	// Accessors and mutators
	conformation::Atom const & atom() const { return base_; }
	void atom( conformation::Atom const & base ) { base_ = base; }

	Vector const & xyz() const { return base_.xyz(); }
	void xyz( Vector const & base ) { base_.xyz(base); }

	Size n_attached_waters() const { return n_attached_waters_; }

	WaterCoords const & waters() const { return waters_; }
	void waters( Size n_attached, WaterCoords const & waters) {
		n_attached_waters_ = n_attached;
		waters_ = waters;
	}

	AtomWeights const & atom_weights() const { return atom_weights_; }
	void atom_weights( AtomWeights const & atom_weights ) { atom_weights_ = atom_weights; }

	/// @brief property required by RotamerTrie class
	inline
	bool is_hydrogen() const { return false; }

	/// @brief send a description of the atom to standard out
	void print() const;

	/// @brief send a description of the atom to an output stream
	void print( std::ostream & os ) const;

	/// @brief compairison operator for sorting
	inline
	bool operator < ( LKBAtom const & other ) const
	{
		// compare the base atom and its coordinate first
		if ( base_.type() != other.base_.type() ) {
			return base_.type() < other.base_.type();
		} else {
			for ( int ii = 0; ii < 3; ++ii ) {
				if ( float(base_.xyz()[ ii ]) != float(other.base_.xyz()[ ii ]) ) {
					return base_.xyz()[ ii ] < other.base_.xyz()[ ii ];
				}
			}
		}

		// Then the coordinates of the attached waters
		if ( n_attached_waters_ != other.n_attached_waters_ ) {
			return n_attached_waters_ < other.n_attached_waters_;
		}

		// because xyzVector does not properly define <
		for ( core::Size i=1; i <= n_attached_waters_; ++i ) {
			for ( core::Size j=0; j<3; ++j ) {
				if ( float(waters_[i][j]) != float(other.waters_[i][j]) ) {
					return (float(waters_[i][j]) < float(other.waters_[i][j]));
				}
			}
		}

		// water coordinates are the same
		for ( core::Size i=1; i <= atom_weights_.size(); ++i ) {
			if ( atom_weights_[ i ] != other.atom_weights_[ i ] ) {
				return atom_weights_[ i ] < other.atom_weights_[ i ];
			}
		}

		// atom weights are the same
		return false;
	}


	/// @brief equality operator for shared-prefix detection -- must be true if
	/// !( a < b ) && !( b < a )
	inline
	bool operator == ( LKBAtom const & other ) const
	{
		if ( base_.type() != other.base_.type() ) {
			return false;
		} else {
			for ( int ii = 0; ii < 3; ++ii ) {
				if ( float(base_.xyz()[ ii ]) != float(other.base_.xyz()[ ii ]) ) {
					return false;
				}
			}
		}

		if ( n_attached_waters_ != other.n_attached_waters_ ) {
			return false;
		}

		// because xyzVector does not properly define <
		for ( core::Size i=1; i <= n_attached_waters_; ++i ) {
			for ( core::Size j=0; j<3; ++j ) {
				if ( float(waters_[i][j]) != float(other.waters_[i][j]) ) {
					return false;
				}
			}
		}

		// waters_ are the same
		for ( core::Size i=1; i <= atom_weights_.size(); ++i ) {
			if ( atom_weights_[ i ] != other.atom_weights_[ i ] ) {
				return false;
			}
		}

		// atom weights are the same

		return true;

	}

private:
	conformation::Atom base_;
	Size n_attached_waters_;
	WaterCoords waters_;
	AtomWeights atom_weights_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

std::ostream & operator << ( std::ostream & os, LKBAtom const & atom );

} // namespace lkbtrie
} // namespace lkball
} // namespace scoring
} // namespace core

#endif
