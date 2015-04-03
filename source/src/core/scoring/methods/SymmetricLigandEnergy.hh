// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/SymmetricLigandEnergy.hh
/// @brief  score for implicit ligand interactions from symmetric geometry
/// @author Will Sheffler


#ifndef INCLUDED_core_scoring_methods_SymmetricLigandEnergy_hh
#define INCLUDED_core_scoring_methods_SymmetricLigandEnergy_hh

// Unit headers
#include <core/scoring/methods/SymmetricLigandEnergy.fwd.hh>


// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>

#include <core/types.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


class SymmetricLigandEnergy : public ContextIndependentOneBodyEnergy {
public:

	typedef ContextIndependentOneBodyEnergy  parent;

	/// @brief ctor
	SymmetricLigandEnergy();

	/// @brief dtor
	virtual ~SymmetricLigandEnergy();

	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////


	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;


	virtual
	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const & sfxn,
		EnergyMap const & emap,
		Vector & F1,
		Vector & F2
	) const;

	/// @brief SymmetricLigandEnergy is context independent; indicates that no
	/// context graphs are required
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const;
virtual
core::Size version() const;


};

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_SymmetricLigandEnergy_HH
