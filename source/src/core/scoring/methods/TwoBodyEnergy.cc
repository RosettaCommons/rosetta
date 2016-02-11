// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/TwoBodyEnergy.cc
/// @brief  Two Body Energy Method base class implementation
/// @author Andrew Leaver-Fay

// Unit Headers
#include <core/scoring/methods/TwoBodyEnergy.hh>

// Package Headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/RotamerSetBase.hh>

#include <basic/Tracer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

static THREAD_LOCAL basic::Tracer TR( "core.scoring.methods.TwoBodyEnergy" );

TwoBodyEnergy::TwoBodyEnergy( EnergyMethodCreatorOP creator ) : parent( creator ) {}

TwoBodyEnergy::~TwoBodyEnergy() {}

bool
TwoBodyEnergy::defines_score_for_residue_pair(
	conformation::Residue const &,
	conformation::Residue const &,
	bool
) const
{
	return true;
}

bool
TwoBodyEnergy::use_extended_residue_pair_energy_interface() const
{
	return false;
}

void
TwoBodyEnergy::residue_pair_energy_ext(
	conformation::Residue const &,
	conformation::Residue const &,
	ResPairMinimizationData const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap &
) const
{
	// APL Noop for now; replace with utility_exit_with_message( "unimplemented extended residue-pair-energy interface method" )
	TR << "WARNING: Unimplemented TwoBodyEnergy::residue_pair_energy_ext()" << std::endl;
}

void
TwoBodyEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const &,
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData &
) const
{
	// Noop
}

/// @details Default return-false implementation.
bool
TwoBodyEnergy::requires_a_setup_for_scoring_for_residue_opportunity( pose::Pose const &  ) const
{
	return false;
}

void
TwoBodyEnergy::setup_for_scoring_for_residue(
	conformation::Residue const & ,
	pose::Pose const & ,
	ScoreFunction const &,
	ResSingleMinimizationData &
) const
{
	// NOOP stub.
	// This should generate an error
	TR << "WARNING: Unimplemented TwoBodyEnergy::setup_for_scoring_for_residue()" << std::endl;
}


bool
TwoBodyEnergy::requires_a_setup_for_derivatives_for_residue_opportunity( pose::Pose const & ) const
{
	return false;
}

void
TwoBodyEnergy::setup_for_derivatives_for_residue(
	conformation::Residue const & ,
	pose::Pose const & ,
	ScoreFunction const &,
	ResSingleMinimizationData &
) const
{
	// NOOP.
	TR << "WARNING: Unimplemented TwoBodyEnergy::setup_for_derivatives_for_residue()" << std::endl;
}

/// @details Default return-false implementation.
bool
TwoBodyEnergy::requires_a_setup_for_scoring_for_residue_pair_opportunity( pose::Pose const &  ) const
{
	return false;
}

void
TwoBodyEnergy::setup_for_scoring_for_residue_pair(
	conformation::Residue const & ,
	conformation::Residue const & ,
	ResSingleMinimizationData const & ,
	ResSingleMinimizationData const & ,
	pose::Pose const & ,
	ScoreFunction const &,
	ResPairMinimizationData &
) const
{
	// NOOP stub.
	// This should generate an error
	TR << "WARNING: Unimplemented TwoBodyEnergy::setup_for_scoring_for_residue_pair()" << std::endl;
}


bool
TwoBodyEnergy::requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose::Pose const & ) const
{
	return false;
}

void
TwoBodyEnergy::setup_for_derivatives_for_residue_pair(
	conformation::Residue const & ,
	conformation::Residue const & ,
	ResSingleMinimizationData const & ,
	ResSingleMinimizationData const & ,
	pose::Pose const & ,
	ScoreFunction const &,
	ResPairMinimizationData &
) const
{
	// NOOP stub.
	TR << "WARNING: Unimplemented TwoBodyEnergy::setup_for_derivatives_for_residue_pair()" << std::endl;
}

void
TwoBodyEnergy::backbone_backbone_energy(
	conformation::Residue const & ,
	conformation::Residue const & ,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap &
) const
{}

void
TwoBodyEnergy::eval_residue_pair_derivatives(
	conformation::Residue const &,
	conformation::Residue const &,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const &,
	pose::Pose const &,
	EnergyMap const &,
	utility::vector1< DerivVectorPair > &,
	utility::vector1< DerivVectorPair > &
) const
{}


void
TwoBodyEnergy::backbone_sidechain_energy(
	conformation::Residue const & ,
	conformation::Residue const & ,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap &
) const
{}


void
TwoBodyEnergy::sidechain_sidechain_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	residue_pair_energy( rsd1, rsd2, pose, sfxn, emap );
}

bool
TwoBodyEnergy::defines_intrares_energy_for_residue(
	conformation::Residue const &
) const
{
	return true;
}

bool
TwoBodyEnergy::use_extended_intrares_energy_interface() const
{
	return false;
}

void
TwoBodyEnergy::eval_intrares_energy_ext(
	conformation::Residue const &,
	ResSingleMinimizationData const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap &
) const
{
	/// APL default noop implementation for this extended interface;
	/// to be replaced with a utility_exit_with_message soon.
	TR << "WARNING: Unimplemented TwoBodyEnergy::eval_intrares_energy_ext()" << std::endl;
}

void
TwoBodyEnergy::setup_for_minimizing_for_residue(
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	ResSingleMinimizationData &
) const
{}

// @brief Low fidelity evaluation of sc/bb + sc/sc energy.
// Do not define in derived class if that class should not be
// used in the packer's bump check phase.
void
TwoBodyEnergy::bump_energy_full(
	conformation::Residue const &,
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap &
) const
{}

void
TwoBodyEnergy::eval_intrares_derivatives(
	conformation::Residue const &,
	ResSingleMinimizationData const &,
	pose::Pose const &,
	EnergyMap const &,
	utility::vector1< DerivVectorPair > &
) const
{
	/// Default noop
}

bool
TwoBodyEnergy::defines_intrares_dof_derivatives( pose::Pose const & ) const
{
	return false;
}

Real
TwoBodyEnergy::eval_intraresidue_dof_derivative(
	conformation::Residue const & ,
	ResSingleMinimizationData const &,
	id::DOF_ID const &,
	id::TorsionID const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap const &
) const
{
	return 0.0;
}

// @brief Low fidelity evaluation of sc/bb energy.
// Do not define in derived class if that class should not be
// used in the packer's bump check phase.
void
TwoBodyEnergy::bump_energy_backbone(
	conformation::Residue const &,
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap &
) const
{}

void
TwoBodyEnergy::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
	using namespace conformation;
	using namespace numeric;

	EnergyMap emap;
	for ( Size ii = 1, ii_end = set1.num_rotamers(); ii <= ii_end; ++ii ) {
		for ( Size jj = 1, jj_end = set2.num_rotamers(); jj <= jj_end; ++jj ) {
			emap.zero();
			residue_pair_energy( *set1.rotamer( ii ), *set2.rotamer( jj ), pose, sfxn, emap );
			energy_table( jj, ii ) += static_cast< core::PackerEnergy > (weights.dot( emap ));
		}
	}
}

void
TwoBodyEnergy::evaluate_rotamer_intrares_energies(
	conformation::RotamerSetBase const & set,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	utility::vector1< core::PackerEnergy > & energies
) const
{
	using namespace conformation;
	using namespace numeric;

	EnergyMap emap;
	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {
		emap.zero();
		eval_intrares_energy( *set.rotamer( ii ), pose, sfxn, emap );
		energies[ ii ] += static_cast< core::PackerEnergy > ( sfxn.weights().dot( emap ) );
	}
}

void
TwoBodyEnergy::evaluate_rotamer_intrares_energy_maps(
	conformation::RotamerSetBase const & set,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	utility::vector1< EnergyMap > & emaps
) const
{
	using namespace conformation;
	using namespace numeric;

	EnergyMap emap;
	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {
		emap.zero();
		eval_intrares_energy( *set.rotamer( ii ), pose, sfxn, emap );
		emaps[ ii ] += emap;
	}
}

void
TwoBodyEnergy::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & residue,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const
{
	EnergyMap emap;
	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {
		emap.zero();
		residue_pair_energy( *set.rotamer( ii ), residue, pose, sfxn, emap );
		energy_vector[ ii ] += static_cast< core::PackerEnergy > (weights.dot( emap ));
	}
}

void
TwoBodyEnergy::evaluate_rotamer_background_energy_maps(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & residue,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & , // weights
	utility::vector1< EnergyMap > & emaps
) const
{
	EnergyMap emap;
	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {
		emap.zero();
		residue_pair_energy( *set.rotamer( ii ), residue, pose, sfxn, emap );
		emaps[ ii ] += emap;
	}
}


}
}
}
