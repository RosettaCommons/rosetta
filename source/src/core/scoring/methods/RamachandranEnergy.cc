// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RamachandranEnergy.cc
/// @brief  Ramachandran energy method class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/methods/RamachandranEnergy.hh>
#include <core/scoring/methods/RamachandranEnergyCreator.hh>

// Package Headers
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/chemical/VariantType.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/id/TorsionID.hh>


// Utility headers
#include <numeric/conversions.hh>
#include <utility/vector1.hh>

//C++ header
#include <stdio.h>

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the RamachandranEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RamachandranEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new RamachandranEnergy;
}

ScoreTypes
RamachandranEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rama );
	return sts;
}


/// ctor
RamachandranEnergy::RamachandranEnergy() :
	parent( new RamachandranEnergyCreator ),
	potential_( ScoringManager::get_instance()->get_Ramachandran() )
{}

/// clone
EnergyMethodOP
RamachandranEnergy::clone() const
{
	return new RamachandranEnergy;
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

///
void
RamachandranEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &pose,
	EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ){
			return;
	}

	if ( rsd.is_protein() &&
				(			(rsd.aa() <= chemical::num_canonical_aas) ||
							(core::chemical::is_canonical_D_aa(rsd.aa()) /*canonical D-amino acids*/) ||
							(rsd.backbone_aa() <= chemical::num_canonical_aas /*noncanonical with canonical template*/)
				)
		) {
			Real rama_score, drama_dphi, drama_dpsi;
			if(potential_.is_normally_connected(rsd)) {
				potential_.eval_rama_score_residue( rsd, rama_score, drama_dphi, drama_dpsi );
			} else {
				potential_.eval_rama_score_residue_nonstandard_connection( pose, rsd, rama_score, drama_dphi, drama_dpsi );
			}
			emap[ rama ] += rama_score;
		}
}

bool
RamachandranEnergy::defines_dof_derivatives( pose::Pose const & ) const
{
	return true;
}

Real
RamachandranEnergy::eval_residue_dof_derivative(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &, // min_data,
	id::DOF_ID const &, //dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const & pose, // pose,
	ScoreFunction const &, // sfxn,
	EnergyMap const & weights
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ){
			return 0.0;
	}

	Real deriv(0.0);
	if ( tor_id.valid() && tor_id.type() == id::BB ) {
		//conformation::Residue const & rsd( pose.residue( tor_id.rsd() ) );
		if (	rsd.is_protein() &&
					(	(rsd.aa() <= chemical::num_canonical_aas) ||
						(rsd.aa()>=core::chemical::aa_dal && rsd.aa()<=core::chemical::aa_dty /*D-amino acids*/) ||
						(rsd.backbone_aa() <= chemical::num_canonical_aas)
					) && tor_id.torsion() <= 2 ) {
			Real rama_score, drama_dphi, drama_dpsi;
			if(potential_.is_normally_connected(rsd)) { //If this residue is connected to the N-1 and N+1 residues
				potential_.eval_rama_score_residue( rsd, rama_score, drama_dphi, drama_dpsi );
				deriv = ( tor_id.torsion() == 1 ? drama_dphi : drama_dpsi );
			} else { //If this residue is connected to things out of sequence
				potential_.eval_rama_score_residue_nonstandard_connection( pose, rsd, rama_score, drama_dphi, drama_dpsi );
				deriv = ( tor_id.torsion() == 1 ? drama_dphi : drama_dpsi );
			}
		}
	}
	// note that the atomtree PHI dofs are in radians
	// use degrees since dE/dangle has angle in denominator
	return numeric::conversions::degrees( weights[ rama ] * deriv );

}



Real
RamachandranEnergy::eval_dof_derivative(
	id::DOF_ID const &,// dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const & pose,
	ScoreFunction const &,// sfxn,
	EnergyMap const & weights
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( pose.residue(tor_id.rsd()).has_variant_type( core::chemical::REPLONLY ) ){
		return 0.0;
	}

	Real deriv(0.0);
	if ( tor_id.valid() && tor_id.type() == id::BB ) {
		conformation::Residue const & rsd( pose.residue( tor_id.rsd() ) );
		if ( rsd.is_protein() &&
					(		(rsd.aa() <= chemical::num_canonical_aas) ||
							(rsd.aa()>=core::chemical::aa_dal && rsd.aa()<=core::chemical::aa_dty /*D-amino acids*/) ||
							(rsd.backbone_aa() <= chemical::num_canonical_aas)
					) && tor_id.torsion() <= 2 ) {
			Real rama_score, drama_dphi, drama_dpsi;
			if(potential_.is_normally_connected(rsd)) { //If this residue is connected to the N-1 and N+1 residues
				potential_.eval_rama_score_residue( rsd, rama_score, drama_dphi, drama_dpsi );
				deriv = ( tor_id.torsion() == 1 ? drama_dphi : drama_dpsi );
			} else { //If this residue is connected to things out of sequence
				potential_.eval_rama_score_residue_nonstandard_connection( pose, rsd, rama_score, drama_dphi, drama_dpsi );
				deriv = ( tor_id.torsion() == 1 ? drama_dphi : drama_dpsi );
			}
		}
	}
	// note that the atomtree PHI dofs are in radians
	// use degrees since dE/dangle has angle in denominator
	return numeric::conversions::degrees( weights[ rama ] * deriv );
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


} // methods
} // scoring
} // core

