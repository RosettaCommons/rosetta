// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/RNA_BaseBaseEnergy.cc
/// @brief  Statistically derived base-base energy
/// @author Rhiju Das


// Unit headers
#include <core/scoring/rna/RNA_BaseBaseEnergy.hh>

// Package headers
#include <core/scoring/rna/RNA_LowResolutionPotential.hh>
#include <core/scoring/rna/RNA_FilteredBaseBaseInfo.hh>
#include <core/scoring/rna/RNA_RawBaseBaseInfo.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyGraph.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>


// Utility headers

#include <ObjexxFCL/format.hh>


// C++


namespace core {
namespace scoring {
namespace rna {

/// c-tor THIS CLASS REQUIRES AN EnergyMethodCreator TO WORK PROPERLY
RNA_BaseBaseEnergy::RNA_BaseBaseEnergy():
	rna_low_resolution_potential_( ScoringManager::get_instance()->get_RNA_LowResolutionPotential() )
{
	/*add_score_type( rna_base_pair );
	add_score_type( rna_base_axis );
	add_score_type( rna_base_stagger );
	add_score_type( rna_base_stack );
	add_score_type( rna_base_stack_axis );*/
}


/// clone
methods::EnergyMethodOP
RNA_BaseBaseEnergy::clone() const
{
	return new RNA_BaseBaseEnergy();
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
///
void
RNA_BaseBaseEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();

	// Following will not repeat any work if any pairwise terms (e.g., base-backbone) have
	// already been calculated... would be better to
	// wrap it inside residue_pair_energy (so it could be used by the packer!).
	// Also maybe this should be unified with RNA_PairwiseLowResolutionEnergy...
	//rna_low_resolution_potential_.update_rna_base_base_interactions( pose );

}


/////////////////////////////////////////////////////////////////////////////
void
RNA_BaseBaseEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
	rna_low_resolution_potential_.update_rna_base_base_interactions( pose );

	rna::RNA_ScoringInfo & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
	rna::RNA_RawBaseBaseInfo const & raw_base_base_info( rna_scoring_info.rna_raw_base_base_info() );
	rna::RNA_FilteredBaseBaseInfo & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );

	rna_filtered_base_base_info.carry_out_filtering( raw_base_base_info );
}


/////////////////////////////////////////////////////////////////////////////
// This all needs to eventually be refactored for the packer.
void
RNA_BaseBaseEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const
{

	//This needs to take care of a bunch of base pair book-keeping.
	// Get Pose cached pairwise "raw" base-base info. This was
	// updated above, in setup_for_scoring
	rna::RNA_ScoringInfo & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
	rna::RNA_RawBaseBaseInfo const & raw_base_base_info( rna_scoring_info.rna_raw_base_base_info() );

	// Create Pose cached non-pairwise, "filtered" base-base info.
	// This forces each base edge to have only one pairing partner.
	rna::RNA_FilteredBaseBaseInfo & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );

	rna_filtered_base_base_info.carry_out_filtering( raw_base_base_info );

	// From filtered base-base info, pull out total scores.
	totals[ rna_base_pair ]        += rna_filtered_base_base_info.get_total_base_pair_score();
	totals[ rna_base_axis ]        += rna_filtered_base_base_info.get_total_base_axis_score();
	totals[ rna_base_stagger ]     += rna_filtered_base_base_info.get_total_base_stagger_score();
	totals[ rna_base_stack ]       += rna_filtered_base_base_info.get_total_base_stack_score();
	totals[ rna_base_stack_axis ]  += rna_filtered_base_base_info.get_total_base_stack_axis_score();



	rna_low_resolution_potential_.finalize( pose ); //Set calculated to false!

	rna_filtered_base_base_info.set_calculated( false );

}

/////////////////////////////////////////////////////////////////////////////
void
RNA_BaseBaseEnergy::eval_atom_derivative(
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
	rna_low_resolution_potential_.eval_atom_derivative_base_base( atom_id, pose, weights, F1, F2 );

	return;
}


/////////////////////////////////////////////////////////////////////////////////////


/// @brief RNA_BaseBaseEnergy distance cutoff
Distance
RNA_BaseBaseEnergy::atomic_interaction_cutoff() const
{
	return 0.0; /// Uh, I don't know.
}

}
}
}
