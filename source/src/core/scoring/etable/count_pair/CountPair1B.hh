// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/count_pair/CountPair1B.hh
/// @brief  Count pair for residue pairs connected with one bond, where the
/// crossover from excluding to counting atom pair interactions is at 3 bonds.
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_scoring_etable_count_pair_CountPair1B_hh
#define INCLUDED_core_scoring_etable_count_pair_CountPair1B_hh


//#include <core/scoring/etable/count_pair/CountPairCrossover3.hh>

#include <core/types.hh>
#include <core/conformation/Atom.hh>

#include <core/conformation/Residue.hh>

#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace etable {
namespace count_pair {

template < class CrossoverBehavior >
class CountPair1B : public CrossoverBehavior
{
public:
public:
	typedef CrossoverBehavior parent;

public:
	CountPair1B(
		conformation::Residue const & res1,
		Size const res1_connect_atom,
		conformation::Residue const & res2,
		Size const res2_connect_atom
	);

	virtual ~CountPair1B() {} // inlined when declared on the stack

	/// @brief function required by templated functions in atom_pair_energy_inline
	inline
	bool
	operator () (
		int const at1,
		int const at2,
		Real & weight,
		Size & path_dist
	) const
	{
		path_dist = 1 + res1_conn_dist_[ at1 ] + res2_conn_dist_[ at2 ];
		return parent::count_at_path_distance( path_dist, weight );
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
	utility::vector1< int > const & res1_conn_dist_;
	utility::vector1< int > const & res2_conn_dist_;

};


/// @brief take a row from the path distances table
/// to retrieve the lower and upper path distances
/// for all atoms in each residue
template < class CrossoverBehavior >
CountPair1B< CrossoverBehavior >::CountPair1B(
	conformation::Residue const & res1,
	Size const res1_connect_atom,
	conformation::Residue const & res2,
	Size const res2_connect_atom
) :
parent(),
res1_conn_dist_( res1.path_distance( res1_connect_atom )),
res2_conn_dist_( res2.path_distance( res2_connect_atom ))
{
}

template < class CrossoverBehavior >
bool
CountPair1B< CrossoverBehavior >::count(
	int const at1,
	int const at2,
	Real & w,
	Size & path_dist
) const
{
	return operator() ( at1, at2, w, path_dist );
}

template < class CrossoverBehavior >
void
CountPair1B< CrossoverBehavior >::residue_atom_pair_energy(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::TableLookupEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy( res1, res2, etable_energy, *this, emap );
}


template < class CrossoverBehavior >
void
CountPair1B< CrossoverBehavior >::residue_atom_pair_energy_sidechain_backbone(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::TableLookupEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_backbone( res1, res2, etable_energy, *this, emap );
}


template < class CrossoverBehavior >
void
CountPair1B< CrossoverBehavior >::residue_atom_pair_energy_sidechain_whole(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::TableLookupEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_whole( res1, res2, etable_energy, *this, emap );
}

template < class CrossoverBehavior >
void
CountPair1B< CrossoverBehavior >::residue_atom_pair_energy_backbone_backbone(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::TableLookupEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_backbone_backbone( res1, res2, etable_energy, *this, emap );
}


template < class CrossoverBehavior >
void
CountPair1B< CrossoverBehavior >::residue_atom_pair_energy_sidechain_sidechain(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::TableLookupEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_sidechain( res1, res2, etable_energy, *this, emap );
}

template < class CrossoverBehavior >
void
CountPair1B< CrossoverBehavior >::residue_atom_pair_energy(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::AnalyticEtableEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy( res1, res2, etable_energy, *this, emap );
}


template < class CrossoverBehavior >
void
CountPair1B< CrossoverBehavior >::residue_atom_pair_energy_sidechain_backbone(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::AnalyticEtableEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_backbone( res1, res2, etable_energy, *this, emap );
}


template < class CrossoverBehavior >
void
CountPair1B< CrossoverBehavior >::residue_atom_pair_energy_sidechain_whole(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::AnalyticEtableEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_whole( res1, res2, etable_energy, *this, emap );
}

template < class CrossoverBehavior >
void
CountPair1B< CrossoverBehavior >::residue_atom_pair_energy_backbone_backbone(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::AnalyticEtableEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_backbone_backbone( res1, res2, etable_energy, *this, emap );
}


template < class CrossoverBehavior >
void
CountPair1B< CrossoverBehavior >::residue_atom_pair_energy_sidechain_sidechain(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::AnalyticEtableEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_sidechain( res1, res2, etable_energy, *this, emap );
}


} // namespace count_pair
} // namespace etable
} // namespace scoring
} // namespace core

#endif
