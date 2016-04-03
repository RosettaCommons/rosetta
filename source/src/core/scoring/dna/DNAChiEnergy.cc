// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dna/DNAChiEnergy.cc
/// @brief  DNA Chi torsion energy method implementation
/// @author Jim Havranek

// Unit headers
#include <core/scoring/dna/DNAChiEnergy.hh>
#include <core/scoring/dna/DNAChiEnergyCreator.hh>
#include <core/scoring/dna/DNABFormPotential.hh>

// Package Headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/scoring/ScoreType.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/id/TorsionID.hh>

// Utility headers
#include <numeric/conversions.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer tr( "core.scoring.dna.DNAChiEnergy" );

namespace core {
namespace scoring {
namespace dna {

methods::EnergyMethodOP
DNAChiEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new DNAChiEnergy );
}

ScoreTypes
DNAChiEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( dna_chi );
	return sts;
}

/// ctor
DNAChiEnergy::DNAChiEnergy() :
	parent( methods::EnergyMethodCreatorOP( new DNAChiEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_DNABFormPotential() )
{}

DNAChiEnergy::~DNAChiEnergy() {}

/// clone
methods::EnergyMethodOP
DNAChiEnergy::clone() const
{
	return methods::EnergyMethodOP( new DNAChiEnergy() );
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////


void
DNAChiEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & ,
	EnergyMap & emap
) const
{
	if ( ! rsd.is_DNA() ) return;

	//  Real this_chi( (rsd.chi())[1] );
	//  tr << "Calculating chi as:  " << this_chi << std::endl;
	Real this_score( 0.0 );
	Real this_deriv( 0.0 );
	potential_.eval_dna_bform_chi_torsion_score_residue( rsd, this_score, this_deriv );
	//  tr << "Calculating dna chi score as:  " << this_score << std::endl;
	emap[ dna_chi ] = this_score;
}


Real
DNAChiEnergy::eval_dof_derivative(
	id::DOF_ID const &,// dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap const & weights
) const
{
	Real this_score( 0.0 );
	Real this_deriv( 0.0 );
	if ( tor_id.valid() && tor_id.type() == id::CHI ) {
		conformation::Residue const & rsd( pose.residue( tor_id.rsd() ) );
		if ( rsd.is_DNA() ) {
			potential_.eval_dna_bform_chi_torsion_score_residue( rsd, this_score, this_deriv );
		}
		//  tr << "Calculating dna chi score as:  " << this_score << std::endl;
	}
	return ( weights[ dna_chi ] * this_deriv );
}

/// @brief DNAChiEnergy is context independent; indicates that no context graphs are required
void
DNAChiEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
) const
{}

core::Size
DNAChiEnergy::version() const
{
	return 1;
}


} // dna
} // scoring
} // core

