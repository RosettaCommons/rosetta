// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/count_pair/CountPairGeneric.hh
/// @brief
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_scoring_etable_count_pair_CountPairGeneric_hh
#define INCLUDED_core_scoring_etable_count_pair_CountPairGeneric_hh

// Unit Headers
#include <core/scoring/etable/count_pair/CountPairGeneric.fwd.hh>

// Package Headers
#include <core/scoring/etable/count_pair/CountPairFunction.hh>

#include <core/types.hh>

#include <core/chemical/ResidueType.fwd.hh>
#include <utility/vector1_bool.hh>


namespace core {
namespace scoring {
namespace etable {
namespace count_pair {

class CountPairGeneric : public CountPairFunction
{
public:
public:
	typedef CountPairFunction parent;

public:
	CountPairGeneric(
		conformation::Residue const & res1,
		conformation::Residue const & res2
	);

	/// @brief Create a count pair object that pretends there exist
	/// a chemical bond between some number of atoms in residue 1 and
	/// some number of atoms in residue2; the bond_pairs vector is
	/// a set of ordered-pairs of atom-indices, where the first
	/// is an atom from restype1 and the second is an atom of restype2.
	CountPairGeneric(
		chemical::ResidueType const & restype1,
		chemical::ResidueType const & restype2,
		utility::vector1< std::pair< Size, Size > > bond_pairs
	);

	void set_crossover( Size );

	virtual ~CountPairGeneric();

	/// @brief function required by templated functions in atom_pair_energy_inline
	inline
	bool
	operator () (
		int const at1,
		int const at2,
		Real & weight,
		Size & minpathdist
	) const
	{
		minpathdist =  path_distance( at1, at2 );

		//std::cout << "CPGeneric: r1= " << r1_.seqpos() << " r2=" << r2_.seqpos();
		//std::cout << " at1: " << r1_.atom_name( at1 ) << " at2: " << r2_.atom_name( at2 );
		//std::cout << " minpathdist= " << minpathdist << " count? ";
		//std::cout << ( minpathdist < crossover_ ? "not counted" : minpathdist == crossover_ ? "half" : "yes")  << std::endl;

		if ( (int) minpathdist < crossover_ ) { return false; }
		if ( (int) minpathdist > crossover_ ) { return true; }
		weight = cp_half;
		return true;
	}

	int
	path_distance( int const at1, int const at2 ) const
	{
		int minpathdist = 10;

		// Check genuine connections
		for ( Size ii = 1; ii <= n_connect_; ++ii ) {
			int ii_pathdist = res1_conn_point_path_dists_[ ii ]->operator[](at1) + res2_conn_point_path_dists_[ ii ]->operator[](at2) + 1;
			if ( ii_pathdist < minpathdist ) minpathdist = ii_pathdist;
		}

		// Check pseudobond connections
		for ( Size ii = 1; ii <= n_pconnect_; ++ii ) {
			int ii_pathdist = res1_pbconn_point_path_dists_[ ii ]->operator[]( at1 )
				+ res2_pbconn_point_path_dists_[ ii ]->operator[]( at2 )
				+ pb_lengths_[ ii ];
			if ( ii_pathdist < minpathdist ) minpathdist = ii_pathdist;
		}
		return minpathdist;
	}

	virtual
	bool
	count(
		int const at1,
		int const at2,
		Real &,
		Size & path_dist
	) const;

	/// Type Resolution Functions ///

	virtual
	void
	residue_atom_pair_energy(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const;


	virtual
	void
	residue_atom_pair_energy_sidechain_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const;


	virtual
	void
	residue_atom_pair_energy_sidechain_whole(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const;


	virtual
	void
	residue_atom_pair_energy_backbone_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const;


	virtual
	void
	residue_atom_pair_energy_sidechain_sidechain(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::TableLookupEvaluator const &,
		EnergyMap &
	) const;

	virtual
	void
	residue_atom_pair_energy(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const;


	virtual
	void
	residue_atom_pair_energy_sidechain_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const;


	virtual
	void
	residue_atom_pair_energy_sidechain_whole(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const;


	virtual
	void
	residue_atom_pair_energy_backbone_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const;


	virtual
	void
	residue_atom_pair_energy_sidechain_sidechain(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::AnalyticEtableEvaluator const &,
		EnergyMap &
	) const;


private:
	//conformation::Residue const & r1_;
	//conformation::Residue const & r2_;

	Size n_connect_;
	utility::vector1< utility::vector1< int > const * > res1_conn_point_path_dists_;
	utility::vector1< utility::vector1< int > const * > res2_conn_point_path_dists_;

	Size n_pconnect_;
	utility::vector1< utility::vector1< int > const * > res1_pbconn_point_path_dists_;
	utility::vector1< utility::vector1< int > const * > res2_pbconn_point_path_dists_;
	utility::vector1< int > pb_lengths_;

	int crossover_;
};

} // namespace count_pair
} // namespace etable
} // namespace scoring
} // namespace core

#endif
