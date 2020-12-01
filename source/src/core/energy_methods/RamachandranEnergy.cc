// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/RamachandranEnergy.cc
/// @brief  Ramachandran energy method class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/energy_methods/RamachandranEnergy.hh>
#include <core/energy_methods/RamachandranEnergyCreator.hh>

// Package Headers
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/chemical/VariantType.hh>

// Project headers
#include <core/id/PartialAtomID.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Residue.functions.hh>
#include <core/pose/Pose.hh>
#include <core/id/TorsionID.hh>


// Utility headers
#include <numeric/conversions.hh>
#include <utility/vector1.hh>

//C++ header
#include <cstdio>
// amw remove me
#include <iomanip>

namespace core {
namespace energy_methods {



/// @details This must return a fresh instance of the RamachandranEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
RamachandranEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< RamachandranEnergy >();
}

core::scoring::ScoreTypes
RamachandranEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( rama );
	return sts;
}


/// ctor
RamachandranEnergy::RamachandranEnergy() :
	parent( utility::pointer::make_shared< RamachandranEnergyCreator >() ),
	potential_( core::scoring::ScoringManager::get_instance()->get_Ramachandran() )
{}

/// clone
core::scoring::methods::EnergyMethodOP
RamachandranEnergy::clone() const
{
	return utility::pointer::make_shared< RamachandranEnergy >();
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////


void
RamachandranEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &pose,
	core::scoring::EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) )  return;

	if ( rsd.is_protein() &&
			(   (rsd.aa() <= chemical::num_canonical_aas) ||
			(core::chemical::is_canonical_D_aa(rsd.aa()) /*canonical D-amino acids*/) ||
			(rsd.backbone_aa() <= chemical::num_canonical_aas /*noncanonical with canonical template*/)
			)
			) {
		Real rama_score, drama_dphi, drama_dpsi;
		if ( potential_.is_normally_connected(rsd) ) {
			potential_.eval_rama_score_residue( rsd, rama_score, drama_dphi, drama_dpsi, (rsd.type().is_achiral_backbone() && rsd.mirrored_relative_to_type()) );
		} else {
			potential_.eval_rama_score_residue_nonstandard_connection( pose, rsd, rama_score, drama_dphi, drama_dpsi, (rsd.type().is_achiral_backbone() && rsd.mirrored_relative_to_type()) );
		}
		//std::cout << "Phi " << std::fixed << std::setprecision( 16 ) << rsd.mainchain_torsions()[1] << std::endl;
		//std::cout << "Psi " << std::fixed << std::setprecision( 16 ) << rsd.mainchain_torsions()[2] << std::endl;
		//std::cout << "Residue energy for " << rsd.name() << " is " << std::fixed << std::setprecision( 16 ) << rama_score << std::endl;
		emap[ core::scoring::rama ] += rama_score;
	}
}

bool
RamachandranEnergy::defines_dof_derivatives( pose::Pose const & ) const
{
	return true;
}

inline
bool
polymeric_termini_incomplete( conformation::Residue res ) {
	for ( Size i = 1; i <= res.n_polymeric_residue_connections(); ++i ) {
		if ( res.connection_incomplete( i ) ) {
			return true;
		}
	}
	return false;
}


/// @details Report which atoms define the score for the %RamachandranEnergy,
/// carefully making sure that if a residue is not elegible for a non-zero
/// score that we don't refer to non-existant atoms (e.g. the atoms before
/// the N-terminus). Use the logic in the Ramachandran class to filter out
/// the termini.
utility::vector1< id::PartialAtomID >
RamachandranEnergy::atoms_with_dof_derivatives( conformation::Residue const & res, pose::Pose const & ) const
{

	utility::vector1< id::PartialAtomID > retlist;

	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( res.has_variant_type( core::chemical::REPLONLY ) )  return retlist;

	//
	if ( ! potential_.defines_score_for_residue( res ) ) {
		return retlist;
	}

	if ( res.is_protein() &&
			( (res.aa() <= chemical::num_canonical_aas) ||
			(res.aa() >= core::chemical::aa_dal && res.aa() <= core::chemical::aa_dty /*D-amino acids*/) ||
			(res.backbone_aa() <= chemical::num_canonical_aas)
			) ) {

		std::set< id::PartialAtomID > atoms;
		for ( Size tor_ind = 1; tor_ind <= 2; ++tor_ind ) {
			conformation::insert_partial_atom_ids_for_mainchain_torsion(
				res, tor_ind, atoms );
		}
		retlist.resize(atoms.size());
		std::copy(atoms.begin(), atoms.end(), retlist.begin());
	}
	return retlist;

}

Real
RamachandranEnergy::eval_residue_dof_derivative(
	conformation::Residue const & rsd,
	core::scoring::ResSingleMinimizationData const &, // min_data,
	id::DOF_ID const &, //dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const & pose, // pose,
	core::scoring::ScoreFunction const &, // sfxn,
	core::scoring::EnergyMap const & weights
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) )  return 0.0;

	// Ignore invalid or non-BB torsions
	if ( !tor_id.valid() || tor_id.type() != id::BB ) return 0.0;

	Real deriv(0.0);
	if ( rsd.is_protein() &&
			( (rsd.aa() <= chemical::num_canonical_aas) ||
			(rsd.aa()>=core::chemical::aa_dal && rsd.aa()<=core::chemical::aa_dty /*D-amino acids*/) ||
			(rsd.backbone_aa() <= chemical::num_canonical_aas)
			) && tor_id.torsion() <= 2
			) {
		Real rama_score, drama_dphi, drama_dpsi;
		if ( potential_.is_normally_connected(rsd) ) { //If this residue is connected to the N-1 and N+1 residues
			potential_.eval_rama_score_residue( rsd, rama_score, drama_dphi, drama_dpsi );
			deriv = ( tor_id.torsion() == 1 ? drama_dphi : drama_dpsi );
		} else { //If this residue is connected to things out of sequence
			potential_.eval_rama_score_residue_nonstandard_connection( pose, rsd, rama_score, drama_dphi, drama_dpsi );
			deriv = ( tor_id.torsion() == 1 ? drama_dphi : drama_dpsi );
		}
		//std::cout << "Phi " << std::fixed << std::setprecision( 16 ) << rsd.mainchain_torsions()[1] << std::endl;
		//std::cout << "Psi " << std::fixed << std::setprecision( 16 ) << rsd.mainchain_torsions()[2] << std::endl;
		//std::cout << "Residue deriv for " << rsd.name() << " is " << deriv << std::endl;
	}

	// note that the atomtree PHI dofs are in radians
	// use degrees since dE/dangle has angle in denominator
	return numeric::conversions::degrees( weights[ core::scoring::rama ] * deriv );

}


Real
RamachandranEnergy::eval_dof_derivative(
	id::DOF_ID const &,// dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &,// sfxn,
	core::scoring::EnergyMap const & weights
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( pose.residue(tor_id.rsd()).has_variant_type( core::chemical::REPLONLY ) ) {
		return 0.0;
	}

	// Ignore invalid or non-BB torsions
	if ( !tor_id.valid() || tor_id.type() != id::BB ) return 0.0;

	Real deriv(0.0);
	conformation::Residue const & rsd( pose.residue( tor_id.rsd() ) );
	if ( rsd.is_protein() &&
			(  (rsd.aa() <= chemical::num_canonical_aas) ||
			(rsd.aa()>=core::chemical::aa_dal && rsd.aa()<=core::chemical::aa_dty /*D-amino acids*/) ||
			(rsd.backbone_aa() <= chemical::num_canonical_aas)
			) && tor_id.torsion() <= 2
			) {
		Real rama_score, drama_dphi, drama_dpsi;
		if ( potential_.is_normally_connected(rsd) ) { //If this residue is connected to the N-1 and N+1 residues
			potential_.eval_rama_score_residue( rsd, rama_score, drama_dphi, drama_dpsi );
			deriv = ( tor_id.torsion() == 1 ? drama_dphi : drama_dpsi );
		} else { //If this residue is connected to things out of sequence
			potential_.eval_rama_score_residue_nonstandard_connection( pose, rsd, rama_score, drama_dphi, drama_dpsi );
			deriv = ( tor_id.torsion() == 1 ? drama_dphi : drama_dpsi );
		}
	}

	// note that the atomtree PHI dofs are in radians
	// use degrees since dE/dangle has angle in denominator
	return numeric::conversions::degrees( weights[ core::scoring::rama ] * deriv );
}

/// @brief Ramachandran Energy is context independent and thus indicates that no context graphs need to
/// be maintained by class Energies
void
RamachandranEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
)
const
{}
core::Size
RamachandranEnergy::version() const
{
	return 1; // Initial versioning
}


} // scoring
} // core

