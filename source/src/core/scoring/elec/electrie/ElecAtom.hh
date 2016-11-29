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

#ifndef INCLUDED_core_scoring_elec_electrie_ElecAtom_hh
#define INCLUDED_core_scoring_elec_electrie_ElecAtom_hh

// Unit Headers
#include <core/scoring/elec/electrie/ElecAtom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.fwd.hh>

// Project Headers
#include <core/types.hh>

// STL Headers
#include <iostream>

// Numceric Headers
#include <numeric/xyzVector.hh>

namespace core {
namespace scoring {
namespace elec {
namespace electrie {

class ElecAtom
{
public:
	ElecAtom();
	ElecAtom( conformation::Residue const & res, Size atom_index );

	~ElecAtom();

	// Accessors and mutators
	conformation::Atom const & atom() const { return base_; }
	void atom( conformation::Atom const & base ) { base_ = base; }

	Vector const & xyz() const { return base_.xyz(); }
	void xyz( Vector const & base ) { base_.xyz(base); }

	utility::vector1< Vector > const & charges() const { return offsite_; }
	void charges(utility::vector1< Vector > const & charges) { offsite_=charges; }

	Real charge() const { return charge_; }
	void charge(Real charge) { charge_=charge; }

	Real frac() const { return frac_; }
	void frac(Real frac) { frac_=frac; }

	void print() const;

	void print( std::ostream & os ) const;

	inline
	bool operator < ( ElecAtom const & other ) const
	{
		if ( offsite_.size() != other.offsite_.size() ) {
			return (offsite_.size() < other.offsite_.size());
		}

		if ( isbb() != other.isbb() ) {
			return (other.isbb());
		}
		if ( charge_ != other.charge_ ) {
			return (charge_ < other.charge_);
		}
		if ( is_hydrogen() != other.is_hydrogen() ) {
			return (other.is_hydrogen());
		}

		// because xyzVector does not properly define <
		for ( core::Size i=1; i<=offsite_.size(); ++i ) {
			for ( core::Size j=0; j<3; ++j ) {
				if ( float(offsite_[i][j]) != float(other.offsite_[i][j]) ) {
					return (float(offsite_[i][j]) < float(other.offsite_[i][j]));
				}
			}
		}

		// offsite_ are the same
		if ( frac_ < other.frac_ ) {
			return true;
		} else if ( frac_ == other.frac_ ) {
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

	inline bool
	isbb() const {
		return isbb_;
	}

	inline void
	isbb(bool isbb)  {
		isbb_=isbb;
	}

	/// @brief property required by RotamerTrie class
	inline
	bool is_hydrogen() const { return is_hydrogen_; }

	/// @brief setter method for data required by RotamerTrie class
	inline
	void is_hydrogen( bool setting ) { is_hydrogen_ = setting; }

	/// @brief equality operator for shared-prefix detection
	inline
	bool operator == ( ElecAtom const & other ) const
	{
		return (
			base_.xyz() == other.base_.xyz() &&
			base_.type() == other.base_.type() &&
			offsite_ == other.offsite_ &&
			isbb_ == other.isbb_ &&
			is_hydrogen_ == other.is_hydrogen_ &&
			charge_ == other.charge_ &&
			frac_ == other.frac_);
	}

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

private:
	conformation::Atom base_;
	utility::vector1< Vector > offsite_;
	Real frac_;
	bool isbb_;
	bool is_hydrogen_;
	Real charge_;
};

std::ostream & operator << ( std::ostream & os, ElecAtom const & atom );

} // namespace lkbtrie
} // namespace lkball
} // namespace scoring
} // namespace core

#endif
