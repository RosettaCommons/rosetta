// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/GenBornEnergy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_methods_GenBornEnergy_hh
#define INCLUDED_core_scoring_methods_GenBornEnergy_hh

// Unit Headers
#include <core/scoring/methods/GenBornEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.hh>
#include <core/scoring/GenBornPotential.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>

#include <utility/vector1.hh>


// Utility headers
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
/**

This is a reimplementation of Jim Havranek's original rosetta++ Gen Born code.
source files: rosetta++/gb_elec*

**/
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


namespace core {
namespace scoring {
namespace methods {


class GenBornEnergy : public ContextDependentLRTwoBodyEnergy  {
public:
	typedef ContextDependentLRTwoBodyEnergy  parent;
public:

	/// for use by ScoringManager
	GenBornEnergy( EnergyMethodOptions const & options );


	GenBornEnergy( GenBornEnergy const & src );


	LongRangeEnergyType long_range_type() const
	{
		return gen_born_lr;
	}

	virtual
	bool
	defines_residue_pair_energy(
		pose::Pose const & pose,
		Size res1,
		Size res2
	) const;

	/// clone
	virtual
	EnergyMethodOP
	clone() const;


	virtual
	void
	setup_for_packing(
		pose::Pose & pose,
		utility::vector1< bool > const & residues_repacking,
		utility::vector1< bool > const &
	) const;


	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;


	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	prepare_rotamers_for_packing(
		pose::Pose const & pose,
		conformation::RotamerSetBase & set
	) const;

	virtual
	void
	update_residue_for_packing(
		pose::Pose &,
		Size resid ) const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	virtual
	void
	evaluate_rotamer_intrares_energies(
		conformation::RotamerSetBase const & set,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		utility::vector1< core::PackerEnergy > & energies
	) const;

	virtual
	void
	evaluate_rotamer_intrares_energy_maps(
		conformation::RotamerSetBase const & set,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		utility::vector1< EnergyMap > & emaps
	) const;

	/// @brief Batch computation of rotamer pair energies.  Need not be overriden in
	/// derived class -- by default, iterates over all pairs of rotamers,
	/// and calls derived class's residue_pair_energy method.  Since short range rotamer pairs
	/// may not need calculation, the default method looks at blocks of residue type pairs
	/// and only calls the residue_pair_energy method if the rotamer pairs are within range
	virtual
	void
	evaluate_rotamer_pair_energies(
		conformation::RotamerSetBase const & set1,
		conformation::RotamerSetBase const & set2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
	) const;


	/// @brief Batch computation of rotamer/background energies.  Need not be overriden
	/// in derived class -- by default, iterates over all rotamers in the set, and calls
	/// derived class's residue_pair_energy method for each one against the background rotamer
	/// Since short range rotamer pairs may not need calculation, the default method
	/// looks at blocks of residue type pairs and only calls the residue_pair_energy method
	/// if the rotamer pairs are within range
	virtual
	void
	evaluate_rotamer_background_energies(
		conformation::RotamerSetBase const & set,
		conformation::Residue const & residue,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		utility::vector1< core::PackerEnergy > & energy_vector
	) const;

	/// @brief Batch computation of rotamer/background energies.  Need not be overriden
	/// in derived class -- by default, iterates over all rotamers in the set, and calls
	/// derived class's residue_pair_energy method for each one against the background rotamer
	/// Since short range rotamer pairs may not need calculation, the default method
	/// looks at blocks of residue type pairs and only calls the residue_pair_energy method
	/// if the rotamer pairs are within range
	virtual
	void
	evaluate_rotamer_background_energy_maps(
		conformation::RotamerSetBase const & set,
		conformation::Residue const & residue,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		utility::vector1< EnergyMap > & emaps
	) const;


	virtual
	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const &,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;

	//  virtual
	//  Distance
	//  atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const;

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	/// this is our own special function
	Real
	packing_interaction_cutoff() const
	{
		return 5.5; // MAGIC NUMBER!!
	}

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	// const-ref to scoring database
	GenBornPotential const & potential_;


	bool const exclude_DNA_DNA_;
	virtual
	core::Size version() const;

};


}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
