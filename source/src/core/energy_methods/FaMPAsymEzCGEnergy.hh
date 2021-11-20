// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/energy_methods/FaMPAsymEzCG.hh
///
/// @brief  Fullatom asymetric EZ potential for CG atoms
/// @details Asymetric EZ potential for CG atoms, from Schramm et al 2012 Structure
///    Last Modified: 7/4/18
///
/// @author  Meghan Franklin (meghanwfranklin@gmail.com)

#ifndef INCLUDED_core_energy_methods_FaMPAsymEzCGEnergy_hh
#define INCLUDED_core_energy_methods_FaMPAsymEzCGEnergy_hh

// Unit Headers
#include <core/energy_methods/FaMPAsymEzCGEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/conformation/Atom.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <ObjexxFCL/FArray2D.hh>
#include <utility/vector1.hh>

namespace core {
namespace energy_methods {

/// @brief Implementation of full-atom Membrane asymmetric EZ potential for the CG atom
///   See Schramm et al 2012 Structure for full details on derivation
///   Statistically derived and calculated derivatives; lookup table based on the Z-coordinate
///   1A bin sizes for -30:30A with a membrane center of 0
class FaMPAsymEzCGEnergy : public core::scoring::methods::ContextIndependentOneBodyEnergy  {

public:

	typedef core::scoring::methods::ContextIndependentOneBodyEnergy  parent;

public:

	// Constructors ////////////////////

	/// @brief Default Constructor
	FaMPAsymEzCGEnergy();

	/// @brief Create a clone of this energy method
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	// Scoring Methods ////////////////

	/// @brief Evaluates the one-body energy for a residue
	void
	residue_energy(
		core::conformation::Residue const & rsd,
		core::pose::Pose const &,
		core::scoring::EnergyMap &
	) const override;

private:

	/// @brief returns the atom name for the atom used to represent the sidechain for
	/// a particular amino acid; this atom was used to derive the statistics this potential
	/// is based on.
	std::string const &
	representative_atom_name( core::chemical::AA const aa ) const;

private:

	ObjexxFCL::FArray2D< core::Real > asymEZ_CG_;

	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_not_required*/ ) const override;

	core::Size version() const override;

};

} // energy_methods
} // core

#endif // INCLUDED_core_energy_methods_FaMPAsymEzCGEnergy_hh
