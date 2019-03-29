// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/AromaticBackboneRestraintEnergy.hh
/// @brief  Aromatic backbone restraint energy method class headers
/// @author Andrew Watkins (amw579@stanford.edu)


#ifndef INCLUDED_core_scoring_methods_AromaticBackboneRestraintEnergy_hh
#define INCLUDED_core_scoring_methods_AromaticBackboneRestraintEnergy_hh

// Unit headers
#include <core/scoring/methods/AromaticBackboneRestraintEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


class AromaticBackboneRestraintEnergy : public ContextIndependentOneBodyEnergy  {
public:
	typedef ContextIndependentOneBodyEnergy  parent;
public:

	/// @brief Constructor.
	///
	AromaticBackboneRestraintEnergy();

	/// @brief Copy constructor.
	///
	AromaticBackboneRestraintEnergy( AromaticBackboneRestraintEnergy const &src );

	/// @brief Clone -- creates a copy and returns an owning pointer to the copy.
	///
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

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }


	/// @brief Evaluate the derivatives for all atoms in this residue.
	///
	virtual
	void
	eval_residue_derivatives(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		pose::Pose const & pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs
	) const;

	/// @brief AromaticBackboneRestraint Energy is context independent and thus indicates that no context graphs need to
	/// be maintained by class Energies
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const;

private:

	virtual
	core::Size version() const;

};

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
