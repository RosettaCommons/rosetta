// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/RNP_LowResStackEnergy.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Kalli Kappel


// Unit headers
#include <core/scoring/rna/RNP_LowResStackEnergy.hh>
#include <core/scoring/rna/RNP_LowResStackData.hh>
#include <core/scoring/rna/RNP_LowResStackEnergyCreator.hh>

// Package headers
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_RawBaseBaseInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <basic/Tracer.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// Utility headers

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

//Auto Headers
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is a statistically derived low-resolution potential for RNA/protein interactions
// For RNA/protein modeling, this is meant to supplement the RNA low-res and protein low-res score
// functions
//
///////////////////////////////////////////////////////////////////////////////////////////////////

static THREAD_LOCAL basic::Tracer TR( "core.scoring.rna.RNP_LowResStackEnergy" );
using namespace core::chemical::rna;
using namespace basic::options;
using namespace basic::options::OptionKeys::score;

namespace core {
namespace scoring {
namespace rna {

/// @details This must return a fresh instance of the RNP_LowResStackEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RNP_LowResStackEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new RNP_LowResStackEnergy );
}

ScoreTypes
RNP_LowResStackEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rnp_stack_xy );
	return sts;
}


/// c-tor
RNP_LowResStackEnergy::RNP_LowResStackEnergy() :
	parent( methods::EnergyMethodCreatorOP( new RNP_LowResStackEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_RNP_LowResStackData() )
{
	//std::cout << "Constructed the RNP stack energy" << std::endl;
}

//clone
methods::EnergyMethodOP
RNP_LowResStackEnergy::clone() const
{
	return methods::EnergyMethodOP( new RNP_LowResStackEnergy );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
//void
//RNP_LowResStackEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & scfxn ) const
//{
//}
//
//void
//RNP_LowResStackEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
//{
//}

//////////////////////////////////////////////////////////////////////////////////////////
void
RNP_LowResStackEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	// Only evaluate these score terms between RNA and protein residues
	if ( !(( rsd1.is_RNA() && rsd2.is_protein() ) || ( rsd1.is_protein() && rsd2.is_RNA() )) ) return;

	Vector rna_centroid;
	Vector protein_centroid;
	numeric::xyzMatrix< core::Real > rna_base_coord_sys;

	if ( rsd1.is_RNA() ) {
		rna_centroid = chemical::rna::get_rna_base_centroid( rsd1 );
		rna_base_coord_sys = chemical::rna::get_rna_base_coordinate_system( rsd1, rna_centroid );
		protein_centroid = rsd2.xyz( "CEN" );
	} else {
		rna_centroid = chemical::rna::get_rna_base_centroid( rsd2 );
		rna_base_coord_sys = chemical::rna::get_rna_base_coordinate_system( rsd2, rna_centroid );
		protein_centroid = rsd1.xyz( "CEN" );
	}
	Vector x_rna = rna_base_coord_sys.col_x();
	Vector y_rna = rna_base_coord_sys.col_y();
	Vector z_rna = rna_base_coord_sys.col_z();

	Vector dist_rna_protein = protein_centroid - rna_centroid;
	Real const dist_x = dot_product( dist_rna_protein, x_rna );
	Real const dist_y = dot_product( dist_rna_protein, y_rna );
	Real const dist_z = dot_product( dist_rna_protein, z_rna );

	// Get the stack score
	Real rnp_stack_xy_score( 0.0 );
	if ( std::abs(dist_z) > 3.0 && std::abs(dist_z) < 6.5 ) {
		potential_.evaluate_rnp_stack_xy_score( rsd1, rsd2, dist_x, dist_y, rnp_stack_xy_score );
	}
	emap[ rnp_stack_xy ] += rnp_stack_xy_score;

}


//////////////////////////////////////////////////////////////////////////////////////////
//void
//RNP_LowResStackEnergy::finalize_total_energy(
// pose::Pose & pose,
// ScoreFunction const &,
// EnergyMap &
//) const {
//
//}

//??
/// @brief RNA_PairwiseLowResolutionEnergy distance cutoff
Distance
RNP_LowResStackEnergy::atomic_interaction_cutoff() const
{
	return 0.0; /// Uh, I don't know.
}

core::Size
RNP_LowResStackEnergy::version() const
{
	return 1; // Initial versioning
}

//etable::count_pair::CountPairFunctionCOP
//RNP_LowResStackEnergy::get_intrares_countpair(
// conformation::Residue const &,
// pose::Pose const &,
// ScoreFunction const &
//) const
//{
// utility_exit_with_message( "FA_ElecEnergy does not define intra - residue pair energies; do not call get_intrares_countpair()" );
// return 0;
//}
//
//etable::count_pair::CountPairFunctionCOP
//RNP_LowResStackEnergy::get_count_pair_function(
// Size const res1,
// Size const res2,
// pose::Pose const & pose,
// ScoreFunction const &
//) const
//{
//}
//
//
//etable::count_pair::CountPairFunctionCOP
//RNP_LowResStackEnergy::get_count_pair_function(
// conformation::Residue const & rsd1,
// conformation::Residue const & rsd2
//) const
//{
//}


} //rna
} //scoring
} //core
