// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/disulfides/CentroidDisulfideEnergy.cc
/// @brief  Centroid Disulfide Energy class implementation
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date   2/4/09


// Unit headers
#include <core/scoring/disulfides/CentroidDisulfideEnergy.hh>
#include <core/scoring/disulfides/CentroidDisulfideEnergyCreator.hh>

// Package headers
#include <core/scoring/disulfides/CentroidDisulfidePotential.hh>
#include <core/scoring/disulfides/CentroidDisulfideEnergyContainer.hh>

// Project headers
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/methods/Methods.hh>
#include <basic/Tracer.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace disulfides {


/// @details This must return a fresh instance of the CentroidDisulfideEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
CentroidDisulfideEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new CentroidDisulfideEnergy( ScoringManager::get_instance()->get_CentroidDisulfidePotential() ) );
}

ScoreTypes
CentroidDisulfideEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( dslfc_cen_dst );
	sts.push_back( dslfc_cb_dst );
	sts.push_back( dslfc_ang );
	sts.push_back( dslfc_cb_dih );
	sts.push_back( dslfc_bb_dih );
	return sts;
}


static thread_local basic::Tracer TR( "core.scoring.disulfides.CentroidDisulfideEnergy" );

CentroidDisulfideEnergy::CentroidDisulfideEnergy(
	CentroidDisulfidePotential const & potential
) :
	parent( methods::EnergyMethodCreatorOP( methods::EnergyMethodCreatorOP( new CentroidDisulfideEnergyCreator ) ) ),
	potential_( potential )
{}

CentroidDisulfideEnergy::~CentroidDisulfideEnergy() {}

// EnergyMethod Methods:

methods::EnergyMethodOP CentroidDisulfideEnergy::clone() const
{
	return methods::EnergyMethodOP( new CentroidDisulfideEnergy( potential_ ) );
}


void CentroidDisulfideEnergy::setup_for_scoring(
		pose::Pose & pose,
		ScoreFunction const & ) const
{
	using namespace methods;

	if ( pose.energies().long_range_container( centroid_disulfide_energy ) == 0 ) {
		CentroidDisulfideEnergyContainerOP dec( new CentroidDisulfideEnergyContainer( pose ) );
		pose.energies().set_long_range_container( centroid_disulfide_energy, dec );
	} else {
		CentroidDisulfideEnergyContainerOP dec = CentroidDisulfideEnergyContainerOP (
				utility::pointer::static_pointer_cast< core::scoring::disulfides::CentroidDisulfideEnergyContainer > ( pose.energies().nonconst_long_range_container( centroid_disulfide_energy ) ));
		dec->update( pose );
	}
}

void CentroidDisulfideEnergy::indicate_required_context_graphs(
		utility::vector1< bool > & ) const
{}


// TwoBodyEnergy Methods:

void CentroidDisulfideEnergy::residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
		) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd1.has_variant_type( core::chemical::REPLONLY ) || rsd2.has_variant_type( core::chemical::REPLONLY ) ){
		return;
	}


	Energy cbcb_distance_score;
	Energy centroid_distance_score;
	Energy cacbcb_angle_1_score;
	Energy cacbcb_angle_2_score;
	Energy cacbcbca_dihedral_score;
	Energy backbone_dihedral_score;

	//Require cysteines
	if ( rsd1.aa() != chemical::aa_cys || rsd2.aa() != chemical::aa_cys ) return;
	//Require Centroid
	if (rsd1.residue_type_set().name() != chemical::CENTROID ||
			rsd2.residue_type_set().name() != chemical::CENTROID )
		return;

	CentroidDisulfideEnergyContainerCOP dec = CentroidDisulfideEnergyContainerCOP (
			utility::pointer::static_pointer_cast< core::scoring::disulfides::CentroidDisulfideEnergyContainer const > ( pose.energies().long_range_container( methods::centroid_disulfide_energy ) ));
	//Require they're bonded
	if ( ! dec->residue_forms_disulfide( rsd1.seqpos() ) ||
			dec->other_neighbor_id( rsd1.seqpos() ) != (Size) rsd2.seqpos() ){
		return;
	}

	potential_.score_disulfide(
			rsd1, rsd2,
			cbcb_distance_score,
			centroid_distance_score,
			cacbcb_angle_1_score,
			cacbcb_angle_2_score,
			cacbcbca_dihedral_score,
			backbone_dihedral_score
			);

	emap[ dslfc_cen_dst ] += centroid_distance_score;
	emap[ dslfc_cb_dst ]  += cbcb_distance_score;
	emap[ dslfc_ang ]     += (cacbcb_angle_1_score + cacbcb_angle_2_score)*.5;
	emap[ dslfc_cb_dih ]  += cacbcbca_dihedral_score;
	emap[ dslfc_bb_dih ]  += backbone_dihedral_score;
}


bool CentroidDisulfideEnergy::defines_intrares_energy( EnergyMap const & ) const
{
	return false;
}


void CentroidDisulfideEnergy::eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
		) const
{}

// LongRangeTwoBodyEnergy methods
methods::LongRangeEnergyType
CentroidDisulfideEnergy::long_range_type() const
{
	return methods::centroid_disulfide_energy;
}


bool CentroidDisulfideEnergy::defines_residue_pair_energy(
		pose::Pose const & pose,
		Size res1,
		Size res2
		) const
{
	using namespace methods;
	if ( ! pose.energies().long_range_container( centroid_disulfide_energy )) return false;

	CentroidDisulfideEnergyContainerCOP dec = CentroidDisulfideEnergyContainerCOP (
			utility::pointer::static_pointer_cast< core::scoring::disulfides::CentroidDisulfideEnergyContainer const > ( pose.energies().long_range_container( centroid_disulfide_energy ) ));
	return dec->disulfide_bonded( res1, res2 );
}
core::Size
CentroidDisulfideEnergy::version() const
{
	return 1; // Initial versioning
}


} // namespace disulfides
} // namespace scoring
} // namespace core

