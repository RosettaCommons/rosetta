// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/PoissonBoltzmannEnergy.hh
/// @brief  Membrane Environment Cbeta Energy
/// @author Bjorn Wallner


#ifndef INCLUDED_core_scoring_methods_PoissonBoltzmannEnergy_HH
#define INCLUDED_core_scoring_methods_PoissonBoltzmannEnergy_HH

// Unit Headers
#include <core/scoring/methods/PoissonBoltzmannEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/PoissonBoltzmannPotential.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <utility/vector1.hh>

#include <string>
#include <map>

namespace core {
namespace scoring {
namespace methods {

class PBLifetimeCache : public basic::datacache::CacheableData{
public:
	PBLifetimeCache();
	virtual ~PBLifetimeCache();
	basic::datacache::CacheableDataOP clone() const;
	void set_charged_residues_map( const std::map<std::string, bool> & charged_residues_map );
	void set_energy_state( const std::string& energy_state );
	// Undefined, commenting out to fix PyRosetta build  void set_pose( const std::string& energy_state, const core::pose::Pose& pose );
	// Undefined, commenting out to fix PyRosetta build  void set_pbp( const std::string& energy_state, core::scoring::PoissonBoltzmannPotentialOP pb );
	void set_conformational_data( const std::string& energy_state,
		const core::pose::Pose& pose,
		core::scoring::PoissonBoltzmannPotentialOP pb );

	bool has_cache( const std::string& energy_state ) const;

	std::map<std::string, bool> & get_charged_residues_map();
	const std::string & get_energy_state() const;
	core::pose::PoseCOP get_pose( const std::string& energy_state );
	core::scoring::PoissonBoltzmannPotentialOP get_pbp( const std::string& energy_state );

private:
	std::map<std::string, bool> charged_residues_map_;
	std::map<std::string, core::pose::PoseCOP> pose_by_state_;
	std::map<std::string, core::scoring::PoissonBoltzmannPotentialOP> pb_by_state_;
	std::string energy_state_;
};


class PoissonBoltzmannEnergy : public ContextIndependentLRTwoBodyEnergy {
public:
	typedef ContextIndependentLRTwoBodyEnergy  parent;
	typedef core::scoring::methods::PBLifetimeCache PBLifetimeCache;
	typedef core::scoring::methods::PBLifetimeCacheOP PBLifetimeCacheOP;
public:


	PoissonBoltzmannEnergy();


	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	methods::LongRangeEnergyType
	long_range_type() const;

	virtual bool defines_intrares_energy( EnergyMap const &  ) const { return true; }

	virtual bool defines_residue_pair_energy(
		pose::Pose const & pose,
		Size res1,
		Size res2
	) const;

	virtual void setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	Real
	revamp_weight_by_burial(
		conformation::Residue const & rsd,
		pose::Pose const & pose
	) const ;
	bool
	residue_in_chains(
		conformation::Residue const & rsd,
		utility::vector1 <Size> chains
	) const ;

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const ;


	virtual void residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	/// virtual
	/// void
	/// residue_energy(
	///  conformation::Residue const & rsd,
	///  pose::Pose const & pose,
	///  EnergyMap & emap
	/// ) const;


	//virtual
	//void
	//finalize_total_energy(
	// pose::Pose & pose,
	// ScoreFunction const &,
	// EnergyMap &
	//) const;

	/// @brief PB Energy is context independent and thus indicates that no context graphs need to
	/// be maintained by class Energies
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:


	/// Compare if two poses are close in fold within tolerance.
	///
	/// To be specific, it returns True if for all a in A and b in B, std::abs(a.xyz - b.xyz) <= tol,
	/// for A and B are sets of all atoms in protein 1 and protein 2, respectively.
	///
	/// @param pose1 A protein's pose
	/// @param pose2 Another protein's pose
	/// @param atom_num The atom number
	/// @param tol   Tolerable distance in Angstrom, >= 0.A
	///
	bool protein_position_equal_within( pose::Pose const & pose1,
		pose::Pose const & pose2,
		Size atom_num,
		Real tolA) const;

	// const-ref to scoring database
	core::Size fixed_residue_;
	core::Real epsilon_;
	mutable scoring::PoissonBoltzmannPotentialOP poisson_boltzmann_potential_;

	// Cached objects for speed
	mutable std::map<std::string, pose::PoseCAP> cached_pose_by_state_;
	mutable std::map<std::string, bool> charged_residues_;
	mutable utility::vector1<Size> charged_chains_;

	virtual
	core::Size version() const;

};


}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
