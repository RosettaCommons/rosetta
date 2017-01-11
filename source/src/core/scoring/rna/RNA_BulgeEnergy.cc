// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/RNA_BulgeEnergy.cc
/// @brief  RNA_Bulge energy method implementation
/// @author Rhiju Das

// Unit headers
#include <core/scoring/rna/RNA_BulgeEnergy.hh>
#include <core/scoring/rna/RNA_BulgeEnergyCreator.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>

// Project headers
#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>


namespace core {
namespace scoring {
namespace rna {

using namespace core::chemical;

/// @details This must return a fresh instance of the RNA_BulgeEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RNA_BulgeEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new RNA_BulgeEnergy );
}

ScoreTypes
RNA_BulgeEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rna_bulge );
	return sts;
}


/// ctor
RNA_BulgeEnergy::RNA_BulgeEnergy() :
	parent( methods::EnergyMethodCreatorOP( new RNA_BulgeEnergyCreator ) ),
	bulge_bonus_( -10.0 ), /*Totally made up for now*/
	rna_bulge_bonus_once_per_loop_( basic::options::option[ basic::options::OptionKeys::score::rna_bulge_bonus_once_per_loop ]() /*default true*/)
{}

RNA_BulgeEnergy::~RNA_BulgeEnergy() {}

/// clone
core::scoring::methods::EnergyMethodOP
RNA_BulgeEnergy::clone() const
{
	return core::scoring::methods::EnergyMethodOP( new RNA_BulgeEnergy );
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////
/// @brief Check for variant types that indicate rsd is an rna_bulge
bool
RNA_BulgeEnergy::is_RNA_bulge( conformation::Residue const & rsd ) const
{
	if ( !rsd.is_RNA() ) return false;
	if ( rsd.has_variant_type( chemical::BULGE ) ) return true;
	if ( rsd.has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) ) return true;
	if ( rsd.has_variant_type( REPLONLY ) ) return false; // I assume this is proper.
	return false;
}


/// @brief Give entropic bonuses for rna_bulge residues
void
RNA_BulgeEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	if ( !is_RNA_bulge( rsd ) ) return;

	if ( rna_bulge_bonus_once_per_loop_ && rsd.seqpos() > 1 ) {
		conformation::Residue const & prev_rsd = pose.residue(rsd.seqpos() - 1);
		if ( rsd.chain() == prev_rsd.chain() && is_RNA_bulge( prev_rsd ) ) {
			return;
		}
	}

	emap[ rna_bulge ] += bulge_bonus_;
}


///////////////////////////////////////////////////////////////////////////////
/// @brief Give entropic bonuses for "missing" residues in the pose
void
RNA_BulgeEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {
	using namespace core::pose::full_model_info;

	// This allows new stepwise monte carlo to encapsulate the physics
	//  modeled in Parin's rna_bulge term. SWM no longer uses VIRTUAL_RNA_RESIDUE or BULGE
	//  variant types, but does have  a full_model_info object that tracks all residues that
	//  need to be built ("missing").
	if ( !full_model_info_defined( pose ) ) return;

	utility::vector1< utility::vector1< Size > > loop_suites;
	Size nmissing = get_number_missing_residues_and_connections( pose, loop_suites );

	if ( rna_bulge_bonus_once_per_loop_ ) {
		// refactor -- have entropic bonus for each loop, but don't increase the
		//  bonus with the number of residues in the loop.
		totals[ rna_bulge ] += bulge_bonus_ * loop_suites.size();
	} else {
		// initial setting -- each missing residue gets a bonus.
		totals[ rna_bulge ] += bulge_bonus_ * nmissing;
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
	return 2; // Adding full_model_info.
}


} //rna
} //scoring
} //core

