// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/methods/AACompositionEnergy.hh
/// @brief Headers for an EnergyMethod that penalizes deviation from a desired amino acid composition.
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.
/// @author Vikram K. Mulligan (vmullig@uw.edu).



#ifndef INCLUDED_core_scoring_methods_AACompositionEnergy_hh
#define INCLUDED_core_scoring_methods_AACompositionEnergy_hh

// Unit headers
#include <core/scoring/methods/AACompositionEnergy.fwd.hh>
#include <core/scoring/methods/AACompositionEnergySetup.fwd.hh>
#include <core/scoring/methods/AACompositionEnergySetup.hh>

// Package headers
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueProperty.hh>

// Project headers
#include <core/types.hh>
#include <map>
#include <string>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

/// @brief AACompositionEnergy, an energy function to penalize stretches of the same residue,
/// derived from base class for EnergyMethods, which are meaningful only on entire structures.
/// These EnergyMethods do all of their work in the "finalize_total_energy" section of score
/// function evaluation.
class AACompositionEnergy : public WholeStructureEnergy {
public:
	typedef WholeStructureEnergy parent;

public:

	/// @brief Default constructor.
	///
	AACompositionEnergy();
	
	/// @brief Copy constructor.
	///
	AACompositionEnergy( AACompositionEnergy const &src );

	/// @brief Default destructor.
	///
	virtual ~AACompositionEnergy();
	
	/// @brief Clone: create a copy of this object, and return an owning pointer
	/// to the copy.
	virtual EnergyMethodOP clone() const;

	/// @brief AACompositionEnergy is context-independent and thus indicates that no context graphs need to be maintained by
	/// class Energies.
	virtual void indicate_required_context_graphs( utility::vector1< bool > &context_graphs_required ) const;

	/// @brief AACompositionEnergy is version 1.0 right now.
	///	
	virtual core::Size version() const;
	
	/// @brief Actually calculate the total energy
	/// @details Called by the scoring machinery.
	virtual void finalize_total_energy( core::pose::Pose & pose, ScoreFunction const &, EnergyMap & totals ) const;
	
	/// @brief Calculate the total energy given a vector of const owning pointers to residues.
	/// @details Called by finalize_total_energy().
	virtual core::Real calculate_aa_composition_energy( utility::vector1< core::conformation::ResidueCOP > const &resvect ) const;
	
	/// @brief Get a summary of all loaded data.
	///
	void report() const;

private:

	/******************
	Private functions:
	******************/

	/// @brief Calculate the total energy based on residue types, given a vector of const owning pointers to residues.
	/// @details Called by calculate_aa_composition_energy().  
	virtual core::Real calculate_energy_by_restype( utility::vector1< core::conformation::ResidueCOP > const &resvect ) const;

	/// @brief Calculate the total energy based on residue properties, given a vector of const owning pointers to residues.
	/// @details Called by calculate_aa_composition_energy().  
	virtual core::Real calculate_energy_by_properties( utility::vector1< core::conformation::ResidueCOP > const &resvect ) const;
	
	/// @brief Get a const-access pointer to the setup helper object.
	///
	inline AACompositionEnergySetupCOP setup_helper_cop() const { return utility::pointer::dynamic_pointer_cast<AACompositionEnergySetup const>(setup_helper_); }
		
	/******************
	Private variables:
	******************/

	/// @brief The helper object that stores all of the data for setting up this scoring function.
	/// @details Initialized on scoreterm initialization.
	AACompositionEnergySetupOP setup_helper_;

};

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
