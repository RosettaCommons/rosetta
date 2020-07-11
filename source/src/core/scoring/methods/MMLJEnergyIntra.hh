// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/MMLJEnergyIntra.hh
/// @brief  molecular mechanics lj energy
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

#ifndef INCLUDED_core_scoring_methods_MMLJEnergyIntra_hh
#define INCLUDED_core_scoring_methods_MMLJEnergyIntra_hh

// Unit headers
#include <core/scoring/methods/MMLJEnergyIntra.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/scoring/NeighborList.hh>
#include <core/scoring/DerivVectorPair.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// C++ headers
#include <iosfwd>

#include <core/scoring/mm/MMLJEnergyTable.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

class NeighborListData : public basic::datacache::CacheableData {
public:

	NeighborListData( scoring::NeighborListOP nblist ) :
		nblist_( nblist )
	{}

	basic::datacache::CacheableDataOP
	clone() const override { return utility::pointer::make_shared< NeighborListData >( *this ); }

	scoring::NeighborListOP
	nblist() const { return nblist_; }

	void
	nblist( scoring::NeighborListOP nblist ) { nblist_ = nblist; }

private:
	//scoring::NeighborList nblist_;
	scoring::NeighborListOP nblist_;
};

typedef utility::pointer::shared_ptr< NeighborListData > NeighborListDataOP;
typedef utility::pointer::shared_ptr< NeighborListData const > NeighborListDataCOP;

class MMLJEnergyIntra : public ContextIndependentTwoBodyEnergy {
public:
	typedef ContextIndependentTwoBodyEnergy  parent;
public:

	/// ctor
	MMLJEnergyIntra();

	/// clone
	EnergyMethodOP
	clone() const override;

	bool
	minimize_in_whole_structure_context( pose::Pose const &  ) const override { return false; }

	void
	setup_for_minimizing_for_residue(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & minmap,
		basic::datacache::BasicDataCache &,
		ResSingleMinimizationData & res_data_cache
	) const override;

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const override;


	void
	eval_intrares_derivatives(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & res_data_cache,
		pose::Pose const & pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs
	) const override;

	bool
	defines_intrares_energy( EnergyMap const & ) const override { return true; }


	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const override;


	bool
	use_extended_intrares_energy_interface() const override { return true; }

	void
	eval_intrares_energy_ext(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & res_data_cache,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const override;

	/// @brief MMLJEnergy does not have an atomic interation threshold
	Distance
	atomic_interaction_cutoff() const override;

	/// @brief MMLJEnergy is context independent; indicates that no context graphs are required
	void
	indicate_required_context_graphs( utility::vector1< bool > & ) const override;

	/// @brief required for neighbor list and to be more lke the ETable
	etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		Size res1,
		Size res2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn
	) const;

	/// @brief required for neighbor list and to be more lke the ETable
	etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn
	) const;

	/// @brief required for neighbor list and to be more lke the ETable
	etable::count_pair::CountPairFunctionOP
	get_intrares_countpair(
		conformation::Residue const & res,
		pose::Pose const & pose,
		ScoreFunction const & sfxn
	) const;

private:
	core::scoring::mm::MMLJEnergyTable const & potential_;
	core::Size version() const override;

};

} // namespace methods
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_methods_MMLJEnergyIntra_HH
