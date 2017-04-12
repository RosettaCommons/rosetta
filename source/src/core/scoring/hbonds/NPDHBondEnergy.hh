// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/NPDHBondEnergy.hh
/// @brief  Non-pairwise-decomposable Hydrogen bond energy method class declaration
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_scoring_hbonds_NPDHBondEnergy_hh
#define INCLUDED_core_scoring_hbonds_NPDHBondEnergy_hh

// Unit Headers
#include <core/scoring/hbonds/NPDHBondEnergy.fwd.hh>

// Package headers
#include <core/scoring/hbonds/hbtrie/HBAtom.hh>
#include <core/scoring/hbonds/hbtrie/HBondTrie.fwd.hh>
#include <core/scoring/hbonds/constants.hh>
//pba

#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/scoring/hbonds/HBondDatabase.fwd.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
#include <utility/vector1.hh>
#include <map>
#include <boost/unordered_map.hpp>

#ifdef PYROSETTA
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#endif

#ifdef    SERIALIZATION
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace hbonds {


/// @brief First pass implementation; somewhat ineffecient.
/// The NPDHbondEnergy will cache a map of all hbonds by donor and acceptor in the
/// Pose during setup_for_scoring, and will then retrieve these values during calls
/// to residue_pair_energy.  Doesn't work properly for packing.
class NPDHBondEnergy : public methods::ContextDependentTwoBodyEnergy  {
public:
	typedef methods::ContextDependentTwoBodyEnergy  parent;
public:

	NPDHBondEnergy( HBondOptions const & opts );
	NPDHBondEnergy( NPDHBondEnergy const & src );

	virtual ~NPDHBondEnergy();

	/// clone
	methods::EnergyMethodOP
	clone() const override;

	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const override;

	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sfxn ) const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const override;

	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const override;

	/// f1 and f2 are zeroed
	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const &,
		ScoreFunction const &,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const override;

	bool minimize_in_whole_structure_context( pose::Pose const & ) const override;

	Distance
	atomic_interaction_cutoff() const override;

	Real
	hydrogen_interaction_cutoff2() const;

	/// @brief HBondEnergy is context sensitive
	void indicate_required_context_graphs(
		utility::vector1< bool > & context_graphs_required ) const override;

	bool
	defines_intrares_energy( EnergyMap const & weights ) const override;

	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const override;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	HBondOptionsCOP options_;
	HBondDatabaseCOP database_;

	core::Size version() const override;

};

} // hbonds
} // scoring
} // core

#ifdef    SERIALIZATION

//CEREAL_FORCE_DYNAMIC_INIT( core_scoring_hbonds_NPDHBondEnergy )

#endif // SERIALIZATION

#endif

