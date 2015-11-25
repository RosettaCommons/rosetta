// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ResidueNeighborList.hh
/// @brief  A container class for use by the Etable and FA_Elec classes for storing
///         lists of atom neighbors
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_ResidueNeighborList_hh
#define INCLUDED_core_scoring_ResidueNeighborList_hh

// Unit headers
#include <core/scoring/ResidueNeighborList.fwd.hh>

// Package headers
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <basic/datacache/CacheableData.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#include <map>

//Auto Headers
namespace core {
namespace scoring {

//enum AtomStatus { uptodate = 0, narrow_ood, wide_ood };

class SmallAtNb
{
public:
	SmallAtNb( Size at1, Size at2, Real weight ) : atomno1_( at1 ), atomno2_( at2 ), weight_( weight ) {}
	Size atomno1() const { return atomno1_; }
	Size atomno2() const { return atomno2_; }
	Real weight() const { return weight_; }
private:
	Size atomno1_;
	Size atomno2_;
	Real weight_;
};

/*class ResidueNblistData : public basic::datacache::CacheableData
{
public:
typedef basic::datacache::CacheableData    parent;
typedef basic::datacache::CacheableDataOP  CacheableDataOP;
typedef basic::datacache::CacheableDataCOP CacheableDataCOP;

typedef utility::in_place_list< AtomNeighbor > AtomNeighbors;
typedef etable::count_pair::CountPairFunctionCOP CountPairFunctionCOP;

public:
ResidueNblistData();
virtual ~ResidueNblistData();

CacheableDataOP clone() const;

void setup_for_intrares_nblist(
CountPairFunctionCOP cpfxn,
Real heavy_heavy_dist_cutoff,
Real heavy_hydrogen_dist_cutoff,
Real hydrogen_hydrogen_dist_cutoff
);

void initialize( conformation::Residue const & res );
void update( conformation::Residue const & res );

int is_H( Size atomno ) const { return is_H_[ atomno ]; }
AtomStatus atom_status( Size atomno ) const { return ood_status_[ atomno ]; }
Size nnarrow_ood() const { return nnarrow_ood_; }
Size nwide_ood() const { return nwide_ood_; }
utility::in_place_list< Size > const & narrow_ood_list() const { return narrow_ood_list_; }
utility::in_place_list< Size > const & wide_ood_list() const { return wide_ood_list_; }
Real narrow_bounding_radius() const { return narrow_bounding_radius_; }
Real wide_bounding_radius() const { return wide_bounding_radius_; }

Vector const & narrow_coord( Size atomno ) const { return narrow_coord_[ atomno ]; }
Vector const & wide_coord( Size atomno ) const { return wide_coord_[ atomno ]; }

AtomNeighbors const & narrow_neighbors( Size atomno ) const { return narrow_[ atomno ]; }
bool intrares_nblist_exists() const { return cpfxn_; }

/// @brief Only use with no-update schemes
utility::vector1< utility::vector1< AtomNeighbor > > const & at_neighbors() const { return at_neighbors_; }

private:
utility::vector1< int >        is_H_;
utility::vector1< AtomStatus > ood_status_;
Size nnarrow_ood_;
Size nwide_ood_;
utility::in_place_list< Size > narrow_ood_list_;
utility::in_place_list< Size > wide_ood_list_;
utility::vector1< Vector >     narrow_coord_;
utility::vector1< Vector >     wide_coord_;

utility::vector1< AtomNeighbors > narrow_;
utility::vector1< AtomNeighbors > wide_;

utility::vector1< utility::vector1< AtomNeighbor > > at_neighbors_; // for no-update schemes

Real narrow_bounding_radius_;
Real wide_bounding_radius_;

CountPairFunctionCOP cpfxn_;
Real ncut2_[ 3 ];
Real wcut2_[ 3 ];
};

/// @brief This class is used to represent for a pair of residues the subset
/// of atom-pair-interactions in a fixed-sequence structure which should be
/// evaluated during score-function evaluation and differentiation.  This
/// class can indeed double as a container for intra-residue atom pair interactions
/// where the residue-1's neighbors for atom i hold the atom's "upper neighbors" --
/// atoms with a larger index -- and residue-2's neighbors for atom i hold
/// it's "lower neighbors".
class ResiduePairNeighborList : public basic::datacache::CacheableData
{
public:
typedef basic::datacache::CacheableData    parent;
typedef basic::datacache::CacheableDataOP  CacheableDataOP;
typedef basic::datacache::CacheableDataCOP CacheableDataCOP;

typedef utility::in_place_list< AtomNeighbor > AtomNeighbors;
typedef etable::count_pair::CountPairFunctionCOP CountPairFunctionCOP;

public:
static const Real narrow_reach; // = { 0.5 };
static const Real wide_reach;   // = { 2.0 };

public:
ResiduePairNeighborList();
virtual ~ResiduePairNeighborList();

CacheableDataOP clone() const;

AtomNeighbors const &
r1_narrow_neighbors( Size atind ) const {
return narrow(0)[ atind ];
}

AtomNeighbors const &
r2_narrow_neighbors( Size atind ) const {
return narrow(1)[ atind ];
}

void initialize_from_residues(
Real heavy_heavy_dist_cutoff,
Real heavy_hydrogen_dist_cutoff,
Real hydrogen_hydrogen_dist_cutoff,
conformation::Residue const & r1,
conformation::Residue const & r2,
ResidueNblistData const & r1dat,
ResidueNblistData const & r2dat,
etable::count_pair::CountPairFunctionCOP cpfxn
);

void update(
Real heavy_heavy_dist_cutoff,
Real heavy_hydrogen_dist_cutoff,
Real hydrogen_hydrogen_dist_cutoff,
conformation::Residue const & r1,
conformation::Residue const & r2,
ResidueNblistData const & r1dat,
ResidueNblistData const & r2dat
);

protected:

AtomNeighbors const &
r1_wide_neighbors( Size atind ) const {
return wide(0)[ atind ];
}

AtomNeighbors const &
r2_wide_neighbors( Size atind ) const {
return wide(1)[ atind ];
}


protected:

/// Assert boundary conditions for index-by-0 arrays here

utility::vector1< AtomNeighbors > &
narrow( int ind ) {
debug_assert( ind == 0 || ind == 1 );
return narrow_[ ind ];
}

utility::vector1< AtomNeighbors > const &
narrow( int ind ) const {
debug_assert( ind == 0 || ind == 1 );
return narrow_[ ind ];
}

utility::vector1< AtomNeighbors > &
wide( int ind ) {
debug_assert( ind == 0 || ind == 1 );
return wide_[ ind ];
}

utility::vector1< AtomNeighbors > const &
wide( int ind ) const {
debug_assert( ind == 0 || ind == 1 );
return wide_[ ind ];
}

public:

/// @brief Only use with no-update schemes
utility::vector1< utility::vector1< AtomNeighbor > > const & r1_at_neighbors() const { return r1_at_neighbors_; }
/// @brief Only use with no-update schemes
utility::vector1< utility::vector1< AtomNeighbor > > const & r2_at_neighbors() const { return r2_at_neighbors_; }

utility::in_place_list< Size > const & r1_nonempty() const { return r1_w_non_empty_narrow_; }

private:
/// @brief debugging functionality: make sure that the atoms in the lists are the ones
/// that should be in the list.
bool
accurate(
Real heavy_heavy_dist_cutoff,
Real heavy_hydrogen_dist_cutoff,
Real hydrogen_hydrogen_dist_cutoff,
conformation::Residue const & r1,
conformation::Residue const & r2,
ResidueNblistData const & r1dat,
ResidueNblistData const & r2dat
) const;


private:
etable::count_pair::CountPairFunctionCOP cpfxn_;

utility::in_place_list< Size > r1_w_non_empty_narrow_;
utility::vector1< AtomNeighbors > narrow_[ 2 ];
utility::vector1< AtomNeighbors > wide_[ 2 ];

utility::vector1< utility::vector1< AtomNeighbor > > r1_at_neighbors_; // for no-update schemes
utility::vector1< utility::vector1< AtomNeighbor > > r2_at_neighbors_; // for no-update schemes

utility::vector1< AtomStatus > update_status_[ 2 ];

Real ncut2_[ 3 ];
Real wcut2_[ 3 ];

};*/

class ResidueNblistData : public basic::datacache::CacheableData
{
public:
	typedef basic::datacache::CacheableData    parent;
	typedef basic::datacache::CacheableDataOP  CacheableDataOP;
	typedef basic::datacache::CacheableDataCOP CacheableDataCOP;

	typedef etable::count_pair::CountPairFunctionCOP CountPairFunctionCOP;

public:
	ResidueNblistData();
	virtual ~ResidueNblistData();

	CacheableDataOP clone() const;

	/// @brief Initialize the residue-nblist; if there are no intra-residue interactions, then provide a null-pointing
	/// count-pair function.
	void initialize(
		conformation::Residue const & res,
		CountPairFunctionCOP cpfxn,
		Real heavy_heavy_dist_cutoff = 0.0,
		Real heavy_hydrogen_dist_cutoff = 0.0,
		Real hydrogen_hydrogen_dist_cutoff = 0.0
	);

	utility::vector1< SmallAtNb > const & atom_neighbors() const { return atom_neighbors_; }

private:
	utility::vector1< SmallAtNb > atom_neighbors_;

};

class ResiduePairNeighborList : public basic::datacache::CacheableData
{
public:
	typedef basic::datacache::CacheableData    parent;
	typedef basic::datacache::CacheableDataOP  CacheableDataOP;
	typedef basic::datacache::CacheableDataCOP CacheableDataCOP;

	//typedef etable::count_pair::CountPairFunctionCOP CountPairFunctionCOP;

public:
	ResiduePairNeighborList();
	virtual ~ResiduePairNeighborList();
	CacheableDataOP clone() const;

	void initialize_from_residues(
		Real vvd2,
		Real hvd2,
		Real hhd2,
		conformation::Residue const & r1,
		conformation::Residue const & r2,
		etable::count_pair::CountPairFunctionCOP cpfxn
	);

	//fpd initialize using count pair representatives
	void initialize_from_residues(
		Real vvd2,
		Real hvd2,
		Real hhd2,
		conformation::Residue const & r1,
		conformation::Residue const & r2,
		etable::count_pair::CountPairFunctionCOP cpfxn,
		std::map<core::Size,core::Size> const &r1_map,
		std::map<core::Size,core::Size> const &r2_map
	);

	utility::vector1< SmallAtNb > const & atom_neighbors() const { return atom_neighbors_; }

private:
	utility::vector1< SmallAtNb > atom_neighbors_;
};

}
}

#endif
