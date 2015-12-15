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
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_hbonds_hbtrie_HBAtom_hh
#define INCLUDED_core_scoring_hbonds_hbtrie_HBAtom_hh

// Unit Headers
#include <core/scoring/hbonds/hbtrie/HBAtom.fwd.hh>

// Package Headers
#include <core/scoring/hbonds/types.hh>

// Project Headers
#include <core/types.hh>

// STL Headers
#include <iostream>

// Numceric Headers
#include <numeric/xyzVector.hh>

namespace core {
namespace scoring {
namespace hbonds {
namespace hbtrie {

class HBAtom
{
public:
	HBAtom();

	~HBAtom();

	// Accessors and mutators
	Vector const & xyz() const { return xyz_; }
	void xyz( Vector const & coord ) { xyz_ = coord; }

	Vector const & base_xyz() const { return base_xyz_; }
	void base_xyz( Vector const & xyz ) { base_xyz_ = xyz; }

	Vector const & base2_xyz() const { return base2_xyz_; }
	void base2_xyz( Vector const & xyz ) { base2_xyz_ = xyz; }

	HBDonChemType hb_don_chem_type() const {debug_assert(  is_hydrogen_ ); return HBDonChemType( hb_chem_type_ ); }
	HBAccChemType hb_acc_chem_type() const {debug_assert( !is_hydrogen_ ); return HBAccChemType( hb_chem_type_ ); }
	void hb_chem_type( int chemtype ) { hb_chem_type_ = chemtype; }

	//int seqpos() const { return seqpos_;}
	//void seqpos( int setting ) { seqpos_ = setting; }


	/// @brief property required by RotamerTrie class
	inline
	bool is_hydrogen() const { return is_hydrogen_; }

	/// @brief setter method for data required by RotamerTrie class
	inline
	void is_hydrogen( bool setting ) { is_hydrogen_ = setting; }

	inline
	bool is_backbone() const { return is_backbone_; }

	inline
	void is_backbone( bool setting ) { is_backbone_ = setting; }

	inline
	bool is_protein() const { return is_protein_; }

	inline
	void is_protein( bool setting ) { is_protein_ = setting; }

	inline
	bool is_dna() const { return is_dna_; }

	inline
	void is_dna( bool setting ) { is_dna_ = setting; }

	/// @brief send a description of the atom to standard out
	void print() const;

	/// @brief send a description of the atom to an output stream
	void print( std::ostream & os ) const;

	/// @brief compairison operator for sorting
	inline
	bool operator < ( HBAtom const & other ) const
	{
		if ( is_hydrogen_ < other.is_hydrogen_ ) {
			return true;
		} else if ( is_hydrogen_ == other.is_hydrogen_ ) {
			if ( hb_chem_type_ < other.hb_chem_type_ ) {
				return true;
			} else if ( hb_chem_type_ == other.hb_chem_type_ ) {
				if ( is_backbone_ < other.is_backbone_ ) {
					return true;
				} else if ( is_backbone_ == other.is_backbone_ ) {
					if ( is_protein_ < other.is_protein_ ) {
						return true;
					} else if ( is_protein_ == other.is_protein_ ) {
						if ( is_dna_ < other.is_dna_ ) {
							return true;
						} else if ( is_dna_ == other.is_dna_ ) {
							//if ( seqpos_ < other.seqpos_ ) {
							// return true;
							//} else if ( seqpos_ == other.seqpos_ ) {
							// Quoting Jack Snoeyink: "Epsilons in code are evil."  But whatcha gonna do?
							// Previously, computed xyz_.distance_squared( other.xyz_ ) and compared against an epsilon
							// In this case, though, you could get points a, b, and c such that a == b, b == c, but a < c.
							// In rare cases this would cause the std::sort() insertion sort step
							// to run off the end of the array/vector and cause a segfault in the trie.
							for ( int ii = 0; ii< 3; ++ii ) {
								if ( float(xyz_[ ii ]) != float(other.xyz_[ ii ]) ) {
									return xyz_[ ii ] < other.xyz_[ ii ];
								}
								if ( float(base_xyz_[ ii ]) != float(other.base_xyz_[ ii ]) ) {
									return base_xyz_[ ii ] < other.base_xyz_[ ii ];
								}

								if ( float(base2_xyz_[ ii ]) != float(other.base2_xyz_[ ii ]) ) {
									return base2_xyz_[ ii ] < other.base2_xyz_[ ii ];
								}
							}
						}
					}
				}
				//}
			}
		}
		return false;
	}


	/// @brief equality operator for shared-prefix detection
	inline
	bool operator == ( HBAtom const & other ) const
	{
		if ( is_hydrogen_ == other.is_hydrogen_ &&
				hb_chem_type_ == other.hb_chem_type_ &&
				is_backbone_ == other.is_backbone_ &&
				is_protein_ == other.is_protein_ &&
				is_dna_ == other.is_dna_
				//&& seqpos_ == other.seqpos_
				) {
			for ( int ii = 0; ii< 3; ++ii ) {
				// Epsilons bad -- see above.
				//Distance diff = std::abs( xyz_[ ii ] - other.xyz_[ ii ] );
				//if ( diff > EPSILON ) {
				if ( float(xyz_[ ii ]) != float(other.xyz_[ ii ]) ) {
					return false;
				}
				//Distance odiff = std::abs( orientation_vector_[ ii ] - other.orientation_vector_[ ii ] );
				//if ( odiff > EPSILON ) {
				if ( float(base_xyz_[ ii ]) != float(other.base_xyz_[ ii ]) ) {
					return false;
				}
				if ( float(base2_xyz_[ ii ]) != float(other.base2_xyz_[ ii ]) ) {
					return false;
				}
			}
		} else {
			return false;
		}
		return true;
	}

	bool non_hbonding_atom() const {
		return hb_chem_type_ == 0;
	}

private:
	Vector xyz_;
	Vector base_xyz_;
	Vector base2_xyz_;

	bool is_hydrogen_;

	bool is_backbone_;

	bool is_protein_;
	bool is_dna_;

	int hb_chem_type_; // an integer either representing an HBDonChemType or an HBAccChemType
	//int seqpos_; // for hbe_classify_BB_by_separation

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

std::ostream & operator << ( std::ostream & os, HBAtom const & atom );

} // namespace hbtrie
} // namespace hbonds
} // namespace scoring
} // namespace core

#endif
