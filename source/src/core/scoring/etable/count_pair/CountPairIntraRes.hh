// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/count_pair/CountPairIntraRes.hh
/// @brief  Count pair for residue pairs connected with one bond, where the
/// crossover from excluding to counting atom pair interactions is at 3 bonds.
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_scoring_etable_count_pair_CountPairIntraRes_hh
#define INCLUDED_core_scoring_etable_count_pair_CountPairIntraRes_hh

// Package Headers
// AUTO-REMOVED #include <core/scoring/etable/count_pair/CountPairFunction.hh>

// Project Headers
#include <core/conformation/Residue.hh>

//#include <core/scoring/etable/count_pair/CountPairCrossover3.hh>

#include <core/types.hh>
#include <core/conformation/Atom.hh>

#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/etable/EtableEnergy.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace etable {
namespace count_pair {

template < class CrossoverBehavior >
class CountPairIntraRes : public CrossoverBehavior
{
public:
	public:
	typedef CrossoverBehavior parent;

public:
	CountPairIntraRes(
		conformation::Residue const & res
	);

	virtual ~CountPairIntraRes();

	///@brief function required by templated functions in atom_pair_energy_inline
	inline
	bool
	operator () (
		int const at1,
		int const at2,
		Real & weight,
		Size & path_dist
	) const
	{
		path_dist = path_dists_[ at1 ][ at2 ];
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
	utility::vector1< utility::vector1< int > > const & path_dists_;

};

/// @brief take a reference to the path distances table
template < class CrossoverBehavior >
CountPairIntraRes< CrossoverBehavior >::CountPairIntraRes(
	conformation::Residue const & res
) :
	parent(),
	path_dists_( res.path_distances() )
{
}

template < class CrossoverBehavior >
CountPairIntraRes< CrossoverBehavior >::~CountPairIntraRes() {}

template < class CrossoverBehavior >
bool
CountPairIntraRes< CrossoverBehavior >::count(
	int const at1,
	int const at2,
	Real & w,
	Size & path_dist
) const
{
	//bool temp = operator() ( at1, at2, w, path_dist );
	//std::cout << "CPIR: " << at1 << "\t" << at2 << "\t" << w << "\t" << path_dist << "\t" << temp << std::endl;
	//return temp;
	return operator() ( at1, at2, w, path_dist );
}


template < class CrossoverBehavior >
void
CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy(
	conformation::Residue const & res,
	conformation::Residue const & ,
	etable::TableLookupEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_intraresidue_atom_pair_energy( res, etable_energy, *this, emap );
}


template < class CrossoverBehavior >
void
CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy_sidechain_backbone(
	conformation::Residue const & ,
	conformation::Residue const & ,
	etable::TableLookupEvaluator const & ,
	EnergyMap &
) const
{
	std::cerr << "Error: illegal call to CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy_sidechain_backbone" << std::endl;
	utility_exit();
}


template < class CrossoverBehavior >
void
CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy_sidechain_whole(
	conformation::Residue const & ,
	conformation::Residue const & ,
	etable::TableLookupEvaluator const & ,
	EnergyMap &
) const
{
	std::cerr << "Error: illegal call to CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy_sidechain_whole" << std::endl;
	utility_exit();
}


template < class CrossoverBehavior >
void
CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy_backbone_backbone(
	conformation::Residue const &,
	conformation::Residue const &,
	etable::TableLookupEvaluator const &,
	EnergyMap &
) const
{
	std::cerr << "Error: illegal call to CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy_backbone_backbone" << std::endl;
	utility_exit();
}


template < class CrossoverBehavior >
void
CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy_sidechain_sidechain(
	conformation::Residue const &,
	conformation::Residue const &,
	etable::TableLookupEvaluator const &,
	EnergyMap &
) const
{
	std::cerr << "Error: illegal call to CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy_sidechain_sidechain" << std::endl;
	utility_exit();
}


template < class CrossoverBehavior >
void
CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy(
	conformation::Residue const & res,
	conformation::Residue const & ,
	etable::AnalyticEtableEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_intraresidue_atom_pair_energy( res, etable_energy, *this, emap );
}


template < class CrossoverBehavior >
void
CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy_sidechain_backbone(
	conformation::Residue const & ,
	conformation::Residue const & ,
	etable::AnalyticEtableEvaluator const & ,
	EnergyMap &
) const
{
	std::cerr << "Error: illegal call to CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy_sidechain_backbone" << std::endl;
	utility_exit();
}


template < class CrossoverBehavior >
void
CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy_sidechain_whole(
	conformation::Residue const & ,
	conformation::Residue const & ,
	etable::AnalyticEtableEvaluator const & ,
	EnergyMap &
) const
{
	std::cerr << "Error: illegal call to CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy_sidechain_whole" << std::endl;
	utility_exit();
}


template < class CrossoverBehavior >
void
CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy_backbone_backbone(
	conformation::Residue const &,
	conformation::Residue const &,
	etable::AnalyticEtableEvaluator const &,
	EnergyMap &
) const
{
	std::cerr << "Error: illegal call to CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy_backbone_backbone" << std::endl;
	utility_exit();
}


template < class CrossoverBehavior >
void
CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy_sidechain_sidechain(
	conformation::Residue const &,
	conformation::Residue const &,
	etable::AnalyticEtableEvaluator const &,
	EnergyMap &
) const
{
	std::cerr << "Error: illegal call to CountPairIntraRes< CrossoverBehavior >::residue_atom_pair_energy_sidechain_sidechain" << std::endl;
	utility_exit();
}


} // namespace count_pair
} // namespace etable
} // namespace scoring
} // namespace core

#endif
