// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/data/RNA_ChemicalMappingEnergy.cc
/// @brief  Statistically derived base-base energy
/// @author Rhiju Das


// Unit headers
#include <core/scoring/rna/data/RNA_ChemicalMappingEnergy.hh>
#include <core/scoring/rna/data/RNA_ChemicalMappingEnergyCreator.hh>
#include <core/scoring/rna/data/RNA_DMS_Potential.hh>

// Package headers
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyGraph.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <basic/Tracer.hh>

///////////////////////////////////////////////////////////////////////////////
//
// Scores based on log-odds of structural features that correlate with
//  chemical reactivities based on 'gold standard' measurements in the Das lab.
//
//
// NOTE: NO DERIVATIVES YET!
//
// -- rhiju, 2014
//
///////////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "core.scoring.rna.data.RNA_ChemicalMappingEnergy" );

namespace core {
namespace scoring {
namespace rna {
namespace data {

methods::EnergyMethodOP
RNA_ChemicalMappingEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {	return new RNA_ChemicalMappingEnergy; }

ScoreTypes
RNA_ChemicalMappingEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rna_chem_map );
	return sts;
}
RNA_ChemicalMappingEnergy::RNA_ChemicalMappingEnergy():
	parent( new RNA_ChemicalMappingEnergyCreator ),
	DMS_potential_( ScoringManager::get_instance()->get_RNA_DMS_Potential() )
{
}

RNA_ChemicalMappingEnergy::~RNA_ChemicalMappingEnergy()
{
}

/// clone
methods::EnergyMethodOP
RNA_ChemicalMappingEnergy::clone() const
{
	return new RNA_ChemicalMappingEnergy();
}

/////////////////////////////////////////////////////////////////////////////
void
RNA_ChemicalMappingEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const
{
	totals[ rna_chem_map ] += calculate_energy( pose );
}
/////////////////////////////////////////////////////////////////////////////
Real
RNA_ChemicalMappingEnergy::calculate_energy( pose::Pose const & pose ) const {

	Real score( 0.0 );

  rna::RNA_ScoringInfo  const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
	rna::data::RNA_DataInfo const & rna_data_info( rna_scoring_info.rna_data_info() );
	RNA_Reactivities const & rna_reactivities = rna_data_info.rna_reactivities();

	if ( rna_reactivities.size() == 0 ) return 0.0;

	DMS_potential_.initialize( pose );
	// get DMS data that needs to be scored.
	// cycle through those data.
	for ( Size n = 1; n <= rna_reactivities.size(); n++ ){
		RNA_Reactivity const & rna_reactivity = rna_reactivities[ n ];
		if ( rna_reactivity.type() == DMS ) score += DMS_potential_.evaluate( pose, rna_reactivity );
	}

	return score;
}

/////////////////////////////////////////////////////////////////////////////
void
RNA_ChemicalMappingEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
 	) const
{

	// updated above, in setup_for_scoring
	return;
}


/////////////////////////////////////////////////////////////////////////////////////


/// @brief RNA_ChemicalMappingEnergy distance cutoff
Distance
RNA_ChemicalMappingEnergy::atomic_interaction_cutoff() const
{
	return 0.0;
}

} //data
} //rna
} //scoring
} //core
