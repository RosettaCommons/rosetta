// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_TorsionEnergy.cc
/// @brief  RNA_Torsion energy method class implementation
/// @author Rhiju Das

// Unit Headers
#include <core/scoring/rna/RNA_TorsionEnergy.hh>
#include <core/scoring/rna/RNA_TorsionEnergyCreator.hh>
#include <core/scoring/rna/RNA_EnergyMethodOptions.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>

// Package Headers
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoringManager.hh>

// Project headers

// Utility headers

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>

// C++

namespace core {
namespace scoring {
namespace rna {


/// @details This must return a fresh instance of the RNA_TorsionEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RNA_TorsionEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new RNA_TorsionEnergy( options.rna_options() ) );
}

ScoreTypes
RNA_TorsionEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rna_torsion );
	sts.push_back( rna_torsion_sc );
	return sts;
}


/// ctor
RNA_TorsionEnergy::RNA_TorsionEnergy( RNA_EnergyMethodOptions const & options,
	RNA_TorsionPotentialOP rna_torsion_potential /* = 0 */ ) :
	parent( methods::EnergyMethodCreatorOP( new RNA_TorsionEnergyCreator ) ),
	options_( options ),
	rna_torsion_potential_( rna_torsion_potential )
{
	if ( rna_torsion_potential_ == 0 ) rna_torsion_potential_ = RNA_TorsionPotentialOP( new RNA_TorsionPotential( options ) );
}

/// clone
methods::EnergyMethodOP
RNA_TorsionEnergy::clone() const
{
	return methods::EnergyMethodOP( new RNA_TorsionEnergy( options_, rna_torsion_potential_ ) );
}


///////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap ) const {

	emap[ rna_torsion ] += rna_torsion_potential_->residue_pair_energy( rsd1, rsd2, pose );

}


///////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap  ) const {

	emap[ rna_torsion ]    += rna_torsion_potential_->eval_intrares_energy( rsd, pose );
	emap[ rna_torsion_sc ] += rna_torsion_potential_->intrares_side_chain_score(); // evaluated at same time as above.

}


///////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{

	rna_torsion_potential_->eval_atom_derivative( id, pose, weights, F1, F2 );

}


void RNA_TorsionEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required */ ) const {}

/// @brief RNA_PairwiseLowResolutionEnergy distance cutoff
Distance
RNA_TorsionEnergy::atomic_interaction_cutoff() const
{
	return 0.0; /// Uh, I don't know.
}

core::Size
RNA_TorsionEnergy::version() const
{
	return 2; // A new torsion potential (integration from Das lab branch -- Aug 2011)
}


} //rna
} //scoring
} //core

