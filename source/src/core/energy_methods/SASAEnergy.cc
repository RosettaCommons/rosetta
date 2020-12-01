// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/SASAEnergy.cc
/// @brief  Power Diagram-derived solvent-accessible surface area energy
/// @author Jim Havranek


// Unit headers
#include <core/energy_methods/SASAEnergy.hh>
#include <core/energy_methods/SASAEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/SASAPotential.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>
// AUTO-REMOVED #include <core/scoring/TenANeighborGraph.hh>
//#include <core/scoring/ContextGraphTypes.hh>

#include <core/scoring/DenseEnergyContainer.hh>

// Project headers
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/conformation/RotamerSetCacheableDataType.hh>

#include <numeric/xyz.io.hh>

#include <utility/vector1.hh>

static basic::Tracer TR( "core.scoring.methods.SASAEnergy" );

// Utility headers


// C++
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


namespace core {
namespace energy_methods {


/// @details This must return a fresh instance of the SASAEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
SASAEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< SASAEnergy >( options );
}

core::scoring::ScoreTypes
SASAEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( fa_sasa );
	return sts;
}


SASAEnergy::SASAEnergy( SASAEnergy const & /*src*/ ) = default;


SASAEnergy::SASAEnergy( core::scoring::methods::EnergyMethodOptions const & options ):
	parent( utility::pointer::make_shared< SASAEnergyCreator >() ),
	potential_( core::scoring::ScoringManager::get_instance()->get_SASAPotential() ),
	exclude_DNA_DNA_( options.exclude_DNA_DNA() )
{}


/// clone
core::scoring::methods::EnergyMethodOP
SASAEnergy::clone() const
{
	return utility::pointer::make_shared< SASAEnergy >( *this );
}

bool
SASAEnergy::defines_residue_pair_energy(
	pose::Pose const &,
	Size,
	Size
) const
{
	return true;
}


///
void
SASAEnergy::setup_for_packing(
	pose::Pose & ,
	utility::vector1< bool > const & , // residues_repacking,
	utility::vector1< bool > const &
) const
{

	// This should utterly fail.  This score term is not for packing.
	TR.Fatal << "You cannot pack with fa_sasa." << std::endl;
	utility_exit_with_message("Attempted to pack with fa_sasa.");
}

void
SASAEnergy::prepare_rotamers_for_packing(
	pose::Pose const & ,
	conformation::RotamerSetBase &
) const
{
	// This should utterly fail.  This score term is not for packing.
	TR.Fatal << "You cannot pack with fa_sasa." << std::endl;
	utility_exit_with_message("Attempted to pack with fa_sasa.");
}


//  void
//  update_residue_for_packing(
//   pose::Pose &,
//   Size /*resid*/ ) const
//  {}
void
SASAEnergy::update_residue_for_packing(
	pose::Pose & ,
	Size
) const
{
}

///
void
SASAEnergy::setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const
{
	core::scoring::methods::LongRangeEnergyType const & lr_type( long_range_type() );

	potential_.setup_for_scoring( pose );

	// create a container
	core::scoring::Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == nullptr ) {
		create_new_lre_container = true;
	} else {
		core::scoring::LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		core::scoring::DenseEnergyContainerOP dec( utility::pointer::static_pointer_cast< core::scoring::DenseEnergyContainer > ( lrc ) );
		if ( dec->size() != pose.size() ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		core::scoring::LREnergyContainerOP new_dec = utility::pointer::make_shared< core::scoring::DenseEnergyContainer >( pose.size(), core::scoring::fa_sasa );
		energies.set_long_range_container( lr_type, new_dec );
	}

}


///
void
SASAEnergy::setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const & ) const
{
	core::scoring::methods::LongRangeEnergyType const & lr_type( long_range_type() );

	potential_.setup_for_scoring( pose );

	// create a container
	core::scoring::Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == nullptr ) {
		create_new_lre_container = true;
	} else {
		core::scoring::LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		core::scoring::DenseEnergyContainerOP dec( utility::pointer::static_pointer_cast< core::scoring::DenseEnergyContainer > ( lrc ) );
		if ( dec->size() != pose.size() ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		core::scoring::LREnergyContainerOP new_dec = utility::pointer::make_shared< core::scoring::DenseEnergyContainer >( pose.size(), core::scoring::fa_sasa );
		energies.set_long_range_container( lr_type, new_dec );
	}

}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

///
void
SASAEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & ,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const
{
	//using core::pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO;
	if ( exclude_DNA_DNA_ && rsd1.is_DNA() && rsd2.is_DNA() ) return;

	emap[ core::scoring::fa_sasa ] += potential_.get_res_res_sasa( rsd1, rsd2 );
}


/////////////////////////////////
// Minimization specific stuff
/////////////////////////////////

void
SASAEnergy::setup_for_minimizing_for_residue(
	conformation::Residue const & ,
	pose::Pose const & , // pose,
	core::scoring::ScoreFunction const &, // scorefxn,
	kinematics::MinimizerMapBase const &, // min_map,
	basic::datacache::BasicDataCache &,
	core::scoring::ResSingleMinimizationData &
) const
{
	return;
}


void
SASAEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const &, // rsd1,
	conformation::Residue const &, // rsd2,
	pose::Pose const &,
	core::scoring::ScoreFunction const &, //scorefxn,
	kinematics::MinimizerMapBase const &, // min_map,
	core::scoring::ResSingleMinimizationData const &,
	core::scoring::ResSingleMinimizationData const &,
	core::scoring::ResPairMinimizationData &
) const
{
	return;
}

bool
SASAEnergy::requires_a_setup_for_scoring_for_residue_opportunity_during_minimization( pose::Pose const & ) const
{
	return false;
}

void
SASAEnergy::setup_for_scoring_for_residue(
	conformation::Residue const & ,
	pose::Pose const & ,
	core::scoring::ScoreFunction const &, // sfxn,
	core::scoring::ResSingleMinimizationData &
) const
{
}

// note: you do not need to define this function or setup_for_derivatives_for_residue because
// the base class already provides both of these no-op implementations; or rather,
// setup_for_derivatives_for_residue won't be called because the base class "requires...opportunity"
// returns false.
bool
SASAEnergy::requires_a_setup_for_derivatives_for_residue_opportunity( pose::Pose const &  ) const
{
	return false;
}


void
SASAEnergy::setup_for_derivatives_for_residue(
	conformation::Residue const & ,
	pose::Pose const & ,
	core::scoring::ScoreFunction const & ,
	core::scoring::ResSingleMinimizationData &,
	basic::datacache::BasicDataCache &
) const
{
}

bool
SASAEnergy::use_extended_residue_pair_energy_interface() const
{
	return true;
}

void
SASAEnergy::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	core::scoring::ResPairMinimizationData const & , // pairdata,
	pose::Pose const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const
{
	if ( exclude_DNA_DNA_ && rsd1.is_DNA() && rsd2.is_DNA() ) return;

	//TR << "Calculating residue pair energy ext" << std::endl;

	emap[ core::scoring::fa_sasa ] += potential_.get_res_res_sasa( rsd1, rsd2 );

}

bool
SASAEnergy::use_extended_intrares_energy_interface() const
{
	return true;
}

void
SASAEnergy::eval_intrares_energy_ext(
	conformation::Residue const & rsd,
	core::scoring::ResSingleMinimizationData const & ,
	pose::Pose const & ,
	core::scoring::ScoreFunction const & ,
	core::scoring::EnergyMap & emap
) const
{
	if ( exclude_DNA_DNA_ && rsd.is_DNA() ) return;

	//TR << "Calculating intraresidue energy ext" << std::endl;

	emap[ core::scoring::fa_sasa ] += potential_.get_res_res_sasa( rsd, rsd );
}


/////////////////////////////////
// End minimization specific stuff
/////////////////////////////////

void
SASAEnergy::evaluate_rotamer_intrares_energies(
	conformation::RotamerSetBase const & ,
	pose::Pose const & ,
	core::scoring::ScoreFunction const & ,
	utility::vector1< core::PackerEnergy > &
) const
{
	return;
}

void
SASAEnergy::evaluate_rotamer_intrares_energy_maps(
	conformation::RotamerSetBase const & ,
	pose::Pose const & ,
	core::scoring::ScoreFunction const & , // sxfn,
	utility::vector1< core::scoring::EnergyMap > &
) const
{
	// using namespace conformation;
	// using namespace numeric;

	return;
}

void
SASAEnergy::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & ,
	conformation::RotamerSetBase const & ,
	pose::Pose const & ,
	core::scoring::ScoreFunction const & , // sfxn,
	core::scoring::EnergyMap const & ,
	ObjexxFCL::FArray2D< core::PackerEnergy > &
) const
{
	return;
}

void
SASAEnergy::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & ,
	conformation::Residue const & ,
	pose::Pose const & ,
	core::scoring::ScoreFunction const & , // sfxn,
	core::scoring::EnergyMap const & ,
	utility::vector1< core::PackerEnergy > &
) const
{

	return;

}

void
SASAEnergy::evaluate_rotamer_background_energy_maps(
	conformation::RotamerSetBase const & ,
	conformation::Residue const & ,
	pose::Pose const & ,
	core::scoring::ScoreFunction const & , // sfxn,
	core::scoring::EnergyMap const & ,
	utility::vector1< core::scoring::EnergyMap > &
) const
{

	return;

}


/// @brief SASAEnergy distance cutoff set to the same cutoff used by EtableEnergy, for now
// Distance
// SASAEnergy::atomic_interaction_cutoff() const
// {
//  return 5.5; /// APL remove this magic number!
// }

/// @brief SASAEnergy requires no context graphs
void
SASAEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{
}

/// @brief SASAEnergy does define intraresidue interactions
bool
SASAEnergy::defines_intrares_energy( core::scoring::EnergyMap const & /*weights*/ ) const
{
	return true;
}

void
SASAEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & ,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const
{

	emap[ core::scoring::fa_sasa ] += potential_.get_res_res_sasa( rsd, rsd  );
}

void
SASAEnergy::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	core::scoring::ResSingleMinimizationData const & ,
	pose::Pose const & pose,
	core::scoring::EnergyMap const & weights,
	utility::vector1< core::scoring::DerivVectorPair > & atom_derivs
) const
{

	Real const factor( weights[ core::scoring::fa_sasa] );

	potential_.eval_residue_pair_derivatives( rsd, rsd, pose, factor,
		atom_derivs, atom_derivs );
}


void
SASAEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	core::scoring::ResSingleMinimizationData const & ,
	core::scoring::ResSingleMinimizationData const & ,
	core::scoring::ResPairMinimizationData const & ,
	pose::Pose const & pose, // provides context
	core::scoring::EnergyMap const & weights,
	utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
	utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs
) const
{

	potential_.eval_residue_pair_derivatives( rsd1, rsd2, pose, weights[ core::scoring::fa_sasa ],
		r1_atom_derivs, r2_atom_derivs );

}


core::Size
SASAEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
