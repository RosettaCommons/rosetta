// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/RNA_PairwiseLowResolutionEnergy.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das


// Unit headers
#include <core/scoring/rna/RNA_PairwiseLowResolutionEnergy.hh>
#include <core/scoring/rna/RNA_PairwiseLowResolutionEnergyCreator.hh>

// Package headers
#include <core/scoring/rna/RNA_LowResolutionPotential.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_RawBaseBaseInfo.hh>
#include <core/scoring/rna/RNA_RawBaseBaseInfo.fwd.hh>
#include <core/scoring/rna/data/RNA_DataInfo.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>


// Utility headers


// C++
using namespace core::chemical::rna;

namespace core {
namespace scoring {
namespace rna {


/// @details This must return a fresh instance of the RNA_PairwiseLowResolutionEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RNA_PairwiseLowResolutionEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new RNA_PairwiseLowResolutionEnergy );
}

ScoreTypes
RNA_PairwiseLowResolutionEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rna_base_pair );
	sts.push_back( rna_base_axis );
	sts.push_back( rna_base_stagger );
	sts.push_back( rna_base_stack );
	sts.push_back( rna_base_stack_axis );
	sts.push_back( rna_base_pair_pairwise );
	sts.push_back( rna_base_axis_pairwise );
	sts.push_back( rna_base_stagger_pairwise );
	sts.push_back( rna_base_stack_pairwise );
	sts.push_back( rna_base_stack_axis_pairwise );
	sts.push_back( rna_base_backbone );
	sts.push_back( rna_backbone_backbone );
	sts.push_back( rna_repulsive );
	sts.push_back( rna_data_base );
	return sts;
}


typedef  numeric::xyzMatrix< Real > Matrix;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Note -- this is not ready for the packer.
// To do:
//  * Have each residue carry with it a virtual centroid, and three virtual atoms
//     that define the base coordinate system (for SPEED)
//  â¢ Actually have the energy calculated in residue_pair_energy rather
//     then the pre-cached values.


/// c-tor
RNA_PairwiseLowResolutionEnergy::RNA_PairwiseLowResolutionEnergy() :
	parent( methods::EnergyMethodCreatorOP( new RNA_PairwiseLowResolutionEnergyCreator ) ),
	rna_low_resolution_potential_( ScoringManager::get_instance()->get_RNA_LowResolutionPotential() )
{}

//clone
methods::EnergyMethodOP
RNA_PairwiseLowResolutionEnergy::clone() const
{
	return methods::EnergyMethodOP( new RNA_PairwiseLowResolutionEnergy );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
RNA_PairwiseLowResolutionEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
	rna_low_resolution_potential_.update_rna_centroid_info( pose );

	rna::RNA_ScoringInfo & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
	rna_raw_base_base_info_ = &rna_scoring_info.rna_raw_base_base_info();

	// Use mini's residue_pair_energy to keep track of the book-keeping, instead of this...
	//rna_low_resolution_potential_.update_rna_base_base_interactions( pose );

	might_be_designing_ = false;
}


void
RNA_PairwiseLowResolutionEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
	rna_low_resolution_potential_.update_rna_centroid_info( pose );

	rna::RNA_ScoringInfo & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
	rna_raw_base_base_info_ = &rna_scoring_info.rna_raw_base_base_info();

	// Just make sure that everything has been calculated. This doesn't use residue_pair_energy for bookkeeping,
	//  it just chugs through all the residue pairs.
	rna_low_resolution_potential_.update_rna_base_base_interactions( pose );

	might_be_designing_ = false;
}

/////////////////////////////////////////////////////////////////
void
RNA_PairwiseLowResolutionEnergy::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const & designing_residues
) const {
	pose.update_residue_neighbors();
	rna_low_resolution_potential_.update_rna_centroid_info( pose );

	rna::RNA_ScoringInfo & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
	rna_raw_base_base_info_ = &rna_scoring_info.rna_raw_base_base_info();

	might_be_designing_ = false;
	for ( Size ii = 1; ii <= designing_residues.size(); ++ii ) {
		if ( designing_residues[ ii ] ) {
			might_be_designing_ = true;
			break;
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////////
//
//
//
void
RNA_PairwiseLowResolutionEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {

	Vector centroid1, centroid2;
	kinematics::Stub stub1, stub2;
	get_centroid_information( rsd1, rsd2, pose, centroid1, centroid2, stub1, stub2 );

	rna_low_resolution_potential_.eval_rna_base_pair_energy( *rna_raw_base_base_info_, rsd1, rsd2, centroid1, centroid2, stub1, stub2 );

	// Note that following are used in packing/designing, although for full de novo modeling
	//  there's a score that filters out base pairings that oversubscribe any one base edge.
	emap[ rna_base_pair_pairwise ]       += rna_base_pair_pairwise_pair_energy( rsd1, rsd2 );
	emap[ rna_base_axis_pairwise ]       += rna_base_axis_pairwise_pair_energy( rsd1, rsd2 );
	emap[ rna_base_stagger_pairwise ]    += rna_base_stagger_pairwise_pair_energy( rsd1, rsd2 );
	emap[ rna_base_stack_pairwise ]      += rna_base_stack_pairwise_pair_energy( rsd1, rsd2 );
	emap[ rna_base_stack_axis_pairwise ] += rna_base_stack_axis_pairwise_pair_energy( rsd1, rsd2 );

	////////////////////////////////////////////////////////////////////////////////
	emap[ rna_base_backbone ]     += rna_low_resolution_potential_.rna_base_backbone_pair_energy( rsd1, rsd2,
		centroid1, centroid2, stub1, stub2 );

	////////////////////////////////////////////////////////////////////////////////
	emap[ rna_backbone_backbone ] += rna_low_resolution_potential_.rna_backbone_backbone_pair_energy( rsd1, rsd2 );
	emap[ rna_repulsive ]         += rna_low_resolution_potential_.rna_repulsive_pair_energy( rsd1, rsd2 );
}

///////////////////////////////////////////////////////////////
void
RNA_PairwiseLowResolutionEnergy::get_centroid_information(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	Vector & centroid1,
	Vector & centroid2,
	kinematics::Stub & stub1,
	kinematics::Stub & stub2
) const {

	// This is not very elegant -- should we consider making the
	// centroid a virtual atom on the residue, and stubs defined by
	// three more virtual atoms?
	rna::RNA_ScoringInfo const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
	rna::RNA_CentroidInfo const & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	if ( might_be_designing_ ) {
		//Recalculate from scratch.
		centroid1 = rna_centroid_info.get_base_centroid( rsd1 );
		centroid2 = rna_centroid_info.get_base_centroid( rsd2 );
		stub1 = rna_centroid_info.get_base_coordinate_system( rsd1, centroid1 );
		stub2 = rna_centroid_info.get_base_coordinate_system( rsd2, centroid2 );
	} else {
		//Use centroid informtion already calculated in setup_for_scoring.
		centroid1 = rna_centroid_info.base_centroids()[ rsd1.seqpos() ];
		centroid2 = rna_centroid_info.base_centroids()[ rsd2.seqpos() ];
		stub1 =  rna_centroid_info.base_stubs()[ rsd1.seqpos() ];
		stub2 =  rna_centroid_info.base_stubs()[ rsd2.seqpos() ];
	}
}


////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
Real
RNA_PairwiseLowResolutionEnergy::rna_base_pair_pairwise_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
) const {
	ObjexxFCL::FArray3D < Real > const & base_pair_array ( rna_raw_base_base_info_->base_pair_array() );
	Real score( 0.0 );
	for ( Size i = 1; i <= NUM_EDGES; i++ ) {
		score += base_pair_array( rsd1.seqpos(), rsd2.seqpos(), i );
		score += base_pair_array( rsd2.seqpos(), rsd1.seqpos(), i );
	}
	return score;
}

/////////////////////////////////////////////////////////////////////////////////
Real
RNA_PairwiseLowResolutionEnergy::rna_base_axis_pairwise_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
) const {
	ObjexxFCL::FArray3D < Real > const & base_axis_array ( rna_raw_base_base_info_->base_axis_array() );
	Real score( 0.0 );
	for ( Size i = 1; i <= NUM_EDGES; i++ ) {
		score += base_axis_array( rsd1.seqpos(), rsd2.seqpos(), i );
		score += base_axis_array( rsd2.seqpos(), rsd1.seqpos(), i );
	}
	return score;
}

/////////////////////////////////////////////////////////////////////////////////
Real
RNA_PairwiseLowResolutionEnergy::rna_base_stagger_pairwise_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
) const {
	ObjexxFCL::FArray3D < Real > const & base_stagger_array ( rna_raw_base_base_info_->base_stagger_array() );
	Real score( 0.0 );
	for ( Size i = 1; i <= NUM_EDGES; i++ ) {
		score += base_stagger_array( rsd1.seqpos(), rsd2.seqpos(), i );
		score += base_stagger_array( rsd2.seqpos(), rsd1.seqpos(), i );
	}
	return score;
}

/////////////////////////////////////////////////////////////////////////////////
Real
RNA_PairwiseLowResolutionEnergy::rna_base_stack_pairwise_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
) const {
	ObjexxFCL::FArray2D < Real > const & stack_array ( rna_raw_base_base_info_->base_stack_array() );
	return stack_array( rsd1.seqpos(), rsd2.seqpos() );
}

/////////////////////////////////////////////////////////////////////////////////
Real
RNA_PairwiseLowResolutionEnergy::rna_base_stack_axis_pairwise_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
) const {
	ObjexxFCL::FArray2D < Real > const & stack_axis_array ( rna_raw_base_base_info_->base_stack_axis_array() );
	return stack_axis_array( rsd1.seqpos(), rsd2.seqpos() );
}


////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_PairwiseLowResolutionEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const & sfxn,
	EnergyMap & totals
) const {

	rna::RNA_ScoringInfo & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );
	rna::RNA_RawBaseBaseInfo & raw_base_base_info( rna_scoring_info.rna_raw_base_base_info() );
	clean_up_rna_two_body_energy_tables( raw_base_base_info, pose );

	//This needs to take care of a bunch of base pair book-keeping.
	// Get Pose cached pairwise "raw" base-base info. This was
	// updated above.
	if ( sfxn.has_nonzero_weight( rna_base_pair ) ||
			sfxn.has_nonzero_weight( rna_base_axis ) ||
			sfxn.has_nonzero_weight( rna_base_stagger ) ||
			sfxn.has_nonzero_weight( rna_base_stack ) ||
			sfxn.has_nonzero_weight( rna_base_stack_axis ) ) {

		// Create Pose cached non-pairwise, "filtered" base-base info.
		// This forces each base edge to have only one pairing partner.
		rna::RNA_FilteredBaseBaseInfo & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );

		// Maybe we should put an if statement.
		rna_filtered_base_base_info.carry_out_filtering( raw_base_base_info );

		// From filtered base-base info, pull out total scores.
		totals[ rna_base_pair ]        = rna_filtered_base_base_info.get_total_base_pair_score();
		totals[ rna_base_axis ]        = rna_filtered_base_base_info.get_total_base_axis_score();
		totals[ rna_base_stagger ]     = rna_filtered_base_base_info.get_total_base_stagger_score();
		totals[ rna_base_stack ]       = rna_filtered_base_base_info.get_total_base_stack_score();
		totals[ rna_base_stack_axis ]  = rna_filtered_base_base_info.get_total_base_stack_axis_score();

		//  rna_filtered_base_base_info.set_calculated( false );
		rna::data::RNA_DataInfo const & rna_data_info( rna_scoring_info.rna_data_info() );
		totals[ rna_data_base ] += rna_filtered_base_base_info.get_data_score( rna_data_info );
	}

	// rna_low_resolution_potential_.finalize( pose );
	// std::cout << "DONE SCORING " << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Hey, if we use the domain_map, can save some computation!!!
void
RNA_PairwiseLowResolutionEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & /*domain_map*/,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	// NOTE -- currently have not put in derivatives for pairwise terms!!!! Would not be too hard actually!
	Vector f1( 0.0 ), f2( 0.0 );

	if ( weights[ rna_base_backbone ] != 0.0 ) {
		rna_low_resolution_potential_.eval_atom_derivative_rna_base_backbone( atom_id, pose, f1, f2 );
		F1 += weights[ rna_base_backbone ] * f1;
		F2 += weights[ rna_base_backbone ] * f2;
	}

	if ( weights[ rna_repulsive ] != 0.0 ) {
		rna_low_resolution_potential_.eval_atom_derivative_rna_repulsive( atom_id, pose, f1, f2 );
		F1 += weights[ rna_repulsive ] * f1;
		F2 += weights[ rna_repulsive ] * f2;
	}

	if ( weights[ rna_backbone_backbone ] != 0.0 ) {
		rna_low_resolution_potential_.eval_atom_derivative_rna_backbone_backbone( atom_id, pose, f1, f2 );
		F1 += weights[ rna_backbone_backbone ] * f1;
		F2 += weights[ rna_backbone_backbone ] * f2;
	}

	// Non-pairwise terms.
	// updated above, in setup_for_scoring
	rna_low_resolution_potential_.eval_atom_derivative_base_base( atom_id, pose, weights, F1, F2 );
}

void
RNA_PairwiseLowResolutionEnergy::clean_up_rna_two_body_energy_tables(
	scoring::rna::RNA_RawBaseBaseInfo & raw_base_base_info,
	pose::Pose & pose
) const {

	//Make sure to zero out any of the base pair, base stack, etc. energies
	// in our special RNA cached energies that are not between neighbors,
	// and were therefore never calculated in residue_pair_energy().
	scoring::rna::RNA_RawBaseBaseInfo raw_base_base_info_save( raw_base_base_info );
	raw_base_base_info.zero();

	scoring::EnergyGraph & energy_graph( pose.energies().energy_graph() );

	for ( Size i = 1; i <= pose.size(); i++ ) {
		for ( utility::graph::Graph::EdgeListIter
				iru  = energy_graph.get_node( i )->upper_edge_list_begin(),
				irue = energy_graph.get_node( i )->upper_edge_list_end();
				iru != irue; ++iru ) {
			Size const j( ( *iru )->get_second_node_ind() );

			raw_base_base_info.copy_values( raw_base_base_info_save, i, j );
			raw_base_base_info.copy_values( raw_base_base_info_save, j, i );
		}
	}
}

/// @brief RNA_PairwiseLowResolutionEnergy distance cutoff
Distance
RNA_PairwiseLowResolutionEnergy::atomic_interaction_cutoff() const
{
	return 0.0; /// Uh, I don't know.
}
core::Size
RNA_PairwiseLowResolutionEnergy::version() const
{
	return 1; // Initial versioning
}

} //rna
} //scoring
} //core
