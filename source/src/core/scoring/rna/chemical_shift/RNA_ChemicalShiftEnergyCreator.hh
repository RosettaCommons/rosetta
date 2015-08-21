// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/scoring/rna/chemical_shift/RNA_ChemicalShiftEnergyCreator.hh
/// @brief  Declaration for the class that connects RNA_ChemicalShiftEnergy with the ScoringManager
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)


#ifndef INCLUDED_core_scoring_rna_chemical_shift_RNA_ChemicalShiftEnergyCreator_HH
#define INCLUDED_core_scoring_rna_chemical_shift_RNA_ChemicalShiftEnergyCreator_HH

#include <core/scoring/methods/EnergyMethodCreator.hh>

namespace core {
namespace scoring {
namespace rna {
namespace chemical_shift {

class RNA_ChemicalShiftEnergyCreator : public methods::EnergyMethodCreator
{
public:
	/// @brief Instantiate a new RNA_ChemicalShiftEnergy
	virtual
	methods::EnergyMethodOP
	create_energy_method(
		methods::EnergyMethodOptions const &
	) const;

	/// @brief Return the set of score types claimed by the EnergyMethod
	/// this EnergyMethodCreator creates in its create_energy_method() function
	virtual
	ScoreTypes
	score_types_for_method() const;

};

} //chemical_shift
} //rna
} //scoring
} //core

#endif
