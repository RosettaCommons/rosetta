// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/P_AA_pp_Energy.cc
/// @brief  Probability of observing an amino acid, given its phi/psi energy method declaration
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/energy_methods/P_AA_pp_Energy.hh>
#include <core/energy_methods/P_AA_pp_EnergyCreator.hh>

// Package headers
#include <core/scoring/P_AA.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/chemical/VariantType.hh>

// Project headers
#include <core/id/PartialAtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Residue.functions.hh>

// Numeric headers
#include <numeric/conversions.hh>

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {



/// @details This must return a fresh instance of the P_AA_pp_Energy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
P_AA_pp_EnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< P_AA_pp_Energy >();
}

core::scoring::ScoreTypes
P_AA_pp_EnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( p_aa_pp );
	return sts;
}


/// ctor
P_AA_pp_Energy::P_AA_pp_Energy() :
	parent( utility::pointer::make_shared< P_AA_pp_EnergyCreator >() ),
	p_aa_( core::scoring::ScoringManager::get_instance()->get_P_AA() )
{}

/// clone
core::scoring::methods::EnergyMethodOP
P_AA_pp_Energy::clone() const
{
	return utility::pointer::make_shared< P_AA_pp_Energy >();
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////


void
P_AA_pp_Energy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	core::scoring::EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) {
		return;
	}

	emap[ core::scoring::p_aa_pp ] += p_aa_.P_AA_pp_energy( rsd );
}


bool
P_AA_pp_Energy::defines_dof_derivatives( pose::Pose const & ) const
{
	return true;
}

utility::vector1< id::PartialAtomID >
P_AA_pp_Energy::atoms_with_dof_derivatives( conformation::Residue const & rsd, pose::Pose const & ) const
{
	if ( ! p_aa_.defines_p_aa_pp_energy_for_res(rsd) ||
			rsd.has_variant_type( core::chemical::REPLONLY ) ) {
		return  utility::vector1< id::PartialAtomID >();
	}

	std::set< id::PartialAtomID > atoms;
	for ( Size tor_ind = 1; tor_ind <= 2; ++tor_ind ) {
		conformation::insert_partial_atom_ids_for_mainchain_torsion(
			rsd, tor_ind, atoms );
	}
	utility::vector1< id::PartialAtomID > retlist(atoms.size());
	std::copy(atoms.begin(), atoms.end(), retlist.begin());
	return retlist;

}

Real
P_AA_pp_Energy::eval_residue_dof_derivative(
	conformation::Residue const & rsd,
	core::scoring::ResSingleMinimizationData const &,// min_data,
	id::DOF_ID const &,// dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const &,// pose,
	core::scoring::ScoreFunction const &,// sfxn,
	core::scoring::EnergyMap const & weights
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) {
		return 0.0;
	}

	if ( ! tor_id.valid() ) return 0.0;
	return numeric::conversions::degrees( weights[ core::scoring::p_aa_pp ] * p_aa_.get_Paa_pp_deriv( rsd, tor_id ));
}


Real
P_AA_pp_Energy::eval_dof_derivative(
	id::DOF_ID const &,// dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &, //sfxn,
	core::scoring::EnergyMap const & weights
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( pose.residue( tor_id.rsd() ).has_variant_type( core::chemical::REPLONLY ) ) {
		return 0.0;
	}

	if ( ! tor_id.valid() ) return 0.0;
	return numeric::conversions::degrees( weights[ core::scoring::p_aa_pp ] * p_aa_.get_Paa_pp_deriv( pose.residue( tor_id.rsd() ), tor_id ));
}

/// @brief P_AA_pp_Energy is context independent; indicates that no
/// context graphs are required
void
P_AA_pp_Energy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{}
core::Size
P_AA_pp_Energy::version() const
{
	return 1; // Initial versioning
}


} // scoring
} // core

