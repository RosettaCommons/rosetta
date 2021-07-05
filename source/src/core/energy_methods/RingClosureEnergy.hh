// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/RingClosureEnergy.hh
/// @brief  Noncanonical ring closure energy method class headers
/// @author Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory


#ifndef INCLUDED_core_energy_methods_RingClosureEnergy_hh
#define INCLUDED_core_energy_methods_RingClosureEnergy_hh

// Unit headers
#include <core/energy_methods/RingClosureEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/MinimizationData.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {



class RingClosureEnergy : public core::scoring::methods::ContextIndependentOneBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentOneBodyEnergy  parent;
public:

	/// @brief Constructor.
	///
	RingClosureEnergy();

	/// @brief Copy constructor.
	///
	RingClosureEnergy( RingClosureEnergy const &src );

	/// @brief Clone -- creates a copy and returns an owning pointer to the copy.
	///
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////


	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const override;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const override { return false; }


	/// @brief Evaluate the derivatives for all atoms in this residue.
	///
	void
	eval_residue_derivatives(
		conformation::Residue const & rsd,
		core::scoring::ResSingleMinimizationData const & min_data,
		pose::Pose const & pose,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & atom_derivs
	) const override;

	/// @brief RingClosure Energy is context independent and thus indicates that no context graphs need to
	/// be maintained by class Energies
	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const override;

private:

	//Private data:

	/// @brief Square of the standard deviation of the harmonic
	/// potential that holds virtual atoms atop real atoms.
	core::Real std_dev_sq_;

	core::Size version() const override;

};

} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
