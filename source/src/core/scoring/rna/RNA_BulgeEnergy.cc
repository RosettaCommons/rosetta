// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_BulgeEnergy.cc
/// @brief  RNA_Bulge energy method implementation
/// @author Rhiju Das

// Unit headers
#include <core/scoring/rna/RNA_BulgeEnergy.hh>
#include <core/scoring/rna/RNA_BulgeEnergyCreator.hh>

// Package Headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>

// Project headers
#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace rna {


/// @details This must return a fresh instance of the RNA_BulgeEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RNA_BulgeEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new RNA_BulgeEnergy;
}

ScoreTypes
RNA_BulgeEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rna_bulge );
	return sts;
}



/// ctor
RNA_BulgeEnergy::RNA_BulgeEnergy() :
	parent( new RNA_BulgeEnergyCreator )
{}

RNA_BulgeEnergy::~RNA_BulgeEnergy() {}

/// clone
core::scoring::methods::EnergyMethodOP
RNA_BulgeEnergy::clone() const
{
	return new RNA_BulgeEnergy;
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

/// @details Allocate the scratch space object on the stack to
/// alieviate thread-safety concerns.  Scratch does not use new.
void
RNA_BulgeEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	EnergyMap & emap
) const
{
	static Real const bulge_bonus = -10.0 /*Totally made up for now*/;

	if ( !rsd.is_RNA() ) return;

	if ( rsd.has_variant_type( chemical::BULGE ) ){
		emap[ rna_bulge ] += bulge_bonus;
	}

	if ( rsd.has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) ){
		emap[ rna_bulge ] += bulge_bonus;
	}

}


/// @brief RNA_BulgeEnergy is context independent; indicates that no context graphs are required
void
RNA_BulgeEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
) const
{}

core::Size
RNA_BulgeEnergy::version() const
{
	return 1; // Initial versioning
}



} //rna
} //scoring
} //core

