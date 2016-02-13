// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/etrie/EtableAtom.hh
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_lkball_lkbtrie_LKBAtom_hh
#define INCLUDED_core_scoring_lkball_lkbtrie_LKBAtom_hh

// Unit Headers
#include <core/scoring/lkball/lkbtrie/LKBAtom.fwd.hh>
#include <core/conformation/Atom.hh>

// Project Headers
#include <core/types.hh>

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

	utility::vector1< Vector > const & waters() const { return waters_; }
	void waters(utility::vector1< Vector > const & waters) { waters_=waters; }

	utility::vector1< Real > const & atom_weights() const { return atom_weights_; }
	void atom_weights(utility::vector1< Real > const & atom_weights) { atom_weights_=atom_weights; }

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
		if ( waters_.size() != other.waters_.size() ) {
			return (waters_.size() < other.waters_.size());
		}

		// because xyzVector does not properly define <
		for ( core::Size i=1; i<=waters_.size(); ++i ) {
			for ( core::Size j=0; j<3; ++j ) {
				if ( float(waters_[i][j]) != float(other.waters_[i][j]) ) {
					return (float(waters_[i][j]) < float(other.waters_[i][j]));
				}
			}
		}

		// waters_ are the same
		if ( atom_weights_ < other.atom_weights_ ) {
			return true;
		} else if ( atom_weights_ == other.atom_weights_ ) {
			if ( base_.type() < other.base_.type() ) {
				return true;
			} else if ( base_.type() == other.base_.type() ) {
				for ( int ii = 0; ii< 3; ++ii ) {
					if ( float(base_.xyz()[ ii ]) != float(other.base_.xyz()[ ii ]) ) {
						return base_.xyz()[ ii ] < other.base_.xyz()[ ii ];
					}
				}
			}
		}

		return false;
	}


	/// @brief equality operator for shared-prefix detection
	inline
	bool operator == ( LKBAtom const & other ) const
	{
		return (
			base_.xyz() == other.base_.xyz() &&
			base_.type() == other.base_.type() &&
			waters_ == other.waters_ &&
			atom_weights_ == other.atom_weights_);
	}

private:
	conformation::Atom base_;
	utility::vector1< Vector > waters_;
	utility::vector1< Real > atom_weights_;
};

std::ostream & operator << ( std::ostream & os, LKBAtom const & atom );

} // namespace lkbtrie
} // namespace lkball
} // namespace scoring
} // namespace core

#endif
