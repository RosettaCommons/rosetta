// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.cc
/// @brief --brief--
/// @author --name-- (--email--)

// Unit headers
#include <--path--/--class--.hh>
#include <--path--/--class--Creator.hh>

#include <core/scoring/methods/TwoBodyEnergy.hh> 
#include <core/scoring/methods/EnergyMethod.hh> 

// Package headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/MinimizationData.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/MinimizerMapBase.hh>

#include <core/scoring/DerivVectorPair.hh>
#include <utility/vector1.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "--namespace_dot--.--class--" );

--namespace--

--class--::--class--():
	core::scoring::methods::TwoBodyEnergy( --class--::class_name() )
{

}

--class--::~--class--(){}

void
--class--::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const;

bool
--class--::defines_score_for_residue_pair(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	bool res_moving_wrt_eachother
) const;

bool
--class--::use_extended_residue_pair_energy_interface() const;

void
--class--::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const;

void
--class--::setup_for_minimizing_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & minmap,
	ResSingleMinimizationData & res_data_cache
) const;

void
--class--::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & minmap,
	ResSingleMinimizationData const & res1_data_cache,
	ResSingleMinimizationData const & res2_data_cache,
	ResPairMinimizationData & data_cache
) const;

bool
--class--::requires_a_setup_for_scoring_for_residue_opportunity( pose::Pose const & pose ) const;

void
--class--::setup_for_scoring_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	ResSingleMinimizationData & min_data
) const;

bool
--class--::requires_a_setup_for_derivatives_for_residue_opportunity( pose::Pose const & pose ) const;

void
--class--::setup_for_derivatives_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	ResSingleMinimizationData & min_data
) const;


bool
--class--::requires_a_setup_for_scoring_for_residue_pair_opportunity( pose::Pose const & pose ) const;

void
--class--::setup_for_scoring_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const & minsingle_data1,
	ResSingleMinimizationData const & minsingle_data2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	ResPairMinimizationData & data_cache
) const;

bool
--class--::requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose::Pose const & pose ) const;

void
--class--::setup_for_derivatives_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const & minsingle_data1,
	ResSingleMinimizationData const & minsingle_data2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	ResPairMinimizationData & data_cache
) const;

void
--class--::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const;

void
--class--::backbone_backbone_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const;

void
--class--::backbone_sidechain_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const;

void
--class--::sidechain_sidechain_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const;

bool
--class--::defines_intrares_energy( EnergyMap const & weights ) const = 0;

void
--class--::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const = 0;

bool
--class--::defines_intrares_energy_for_residue(
	conformation::Residue const & res
) const;

bool
--class--::use_extended_intrares_energy_interface() const;

void
--class--::eval_intrares_energy_ext(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & data_cache,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const;

void
--class--::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & min_data,
	pose::Pose const & pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const;

bool
--class--::defines_intrares_dof_derivatives( pose::Pose const & p ) const;

Real
--class--::eval_intraresidue_dof_derivative(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & min_data,
	id::DOF_ID const & dof_id,
	id::TorsionID const & torsion_id,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights
) const;

void
--class--::bump_energy_full(
	conformation::Residue const &,
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap &
) const;

virtual
void
--class--::bump_energy_backbone(
	conformation::Residue const &,
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap &
) const;

void
--class--::evaluate_rotamer_intrares_energies(
	conformation::RotamerSetBase const & set,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	utility::vector1< core::PackerEnergy > & energies
) const;

void
--class--::evaluate_rotamer_intrares_energy_maps(
	conformation::RotamerSetBase const & set,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	utility::vector1< EnergyMap > & emaps
) const;

void
--class--::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const;

void
--class--::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & residue,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const;

void
--class--::evaluate_rotamer_background_energy_maps(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & residue,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	utility::vector1< EnergyMap > & emaps
) const;

--end_namespace--

