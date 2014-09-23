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
#include <core/scoring/rna/data/RNA_DMS_LowResolutionPotential.hh>

// Package headers
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>

#include <basic/Tracer.hh>

///////////////////////////////////////////////////////////////////////////////
//
// Scores based on log-odds of structural features that correlate with
//  chemical reactivities based on 'gold standard' measurements in the Das lab.
//
// NOTE: NO DERIVATIVES YET!
//
// -- rhiju, 2014
//
///////////////////////////////////////////////////////////////////////////////

static thread_local basic::Tracer TR( "core.scoring.rna.data.RNA_ChemicalMappingEnergy" );

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
	sts.push_back( rna_chem_map_lores );
	return sts;
}
RNA_ChemicalMappingEnergy::RNA_ChemicalMappingEnergy():
	parent( methods::EnergyMethodCreatorOP( new RNA_ChemicalMappingEnergyCreator ) ),
	DMS_potential_( ScoringManager::get_instance()->get_RNA_DMS_Potential() ),
	DMS_low_resolution_potential_( ScoringManager::get_instance()->get_RNA_DMS_LowResolutionPotential() )
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
	ScoreFunction const & scorefxn,
	EnergyMap & totals
) const
{
	bool const rna_base_pair_computed = scorefxn.has_nonzero_weight( rna_base_pair ); // assumes this was *earlier* in finalize_total_energy loop.
	if ( scorefxn.has_nonzero_weight( rna_chem_map ) )       totals[ rna_chem_map ]       += calculate_energy( pose, false /*use_low_res*/, rna_base_pair_computed );
	if ( scorefxn.has_nonzero_weight( rna_chem_map_lores ) ) totals[ rna_chem_map_lores ] += calculate_energy( pose, true  /*use_low_res*/, rna_base_pair_computed );
}
/////////////////////////////////////////////////////////////////////////////
Real
RNA_ChemicalMappingEnergy::calculate_energy( pose::Pose & pose,
																						 bool const use_low_res /* = false */,
																						 bool const rna_base_pair_computed /* = false */ ) const {

	Real score( 0.0 );

	pose::full_model_info::make_sure_full_model_info_is_setup( pose );
  rna::RNA_ScoringInfo  const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
	rna::data::RNA_DataInfo const & rna_data_info( rna_scoring_info.rna_data_info() );
	RNA_Reactivities const & rna_reactivities = rna_data_info.rna_reactivities();

	if ( rna_reactivities.size() == 0 ) return 0.0;

	if ( use_low_res ) {
		DMS_low_resolution_potential_.initialize( pose, rna_base_pair_computed );
	} else {
		DMS_potential_.initialize( pose );
	}

	// get DMS data that needs to be scored.
	// cycle through those data.
	for ( Size n = 1; n <= rna_reactivities.size(); n++ ){
		RNA_Reactivity const & rna_reactivity = rna_reactivities[ n ];
		if ( rna_reactivity.type() == DMS  ) 	{
			if ( use_low_res ) {
				score += DMS_low_resolution_potential_.evaluate( pose, rna_reactivity );
			}	else {
				score += DMS_potential_.evaluate( pose, rna_reactivity );
			}
		}
	}

	return score;
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
