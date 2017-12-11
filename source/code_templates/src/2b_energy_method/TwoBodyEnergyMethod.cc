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

static basic::Tracer TR( "--namespace_dot--.--class--" );

--namespace--

using namespace core::scoring;

--class--::--class--():
	core::scoring::methods::TwoBodyEnergy( --class--::class_name() )
{

}

--class--::~--class--(){}

void
--class--::residue_pair_energy(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const;

bool
--class--::defines_score_for_residue_pair(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2,
	bool res_moving_wrt_eachother
) const;

bool
--class--::use_extended_residue_pair_energy_interface() const;

void
--class--::residue_pair_energy_ext(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	ResPairMinimizationData const & min_data,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const;

void
--class--::setup_for_minimizing_for_residue(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	core::kinematics::MinimizerMapBase const & minmap,
	ResSingleMinimizationData & res_data_cache
) const;

void
--class--::setup_for_minimizing_for_residue_pair(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	core::kinematics::MinimizerMapBase const & minmap,
	ResSingleMinimizationData const & res1_data_cache,
	ResSingleMinimizationData const & res2_data_cache,
	ResPairMinimizationData & data_cache
) const;

bool
--class--::requires_a_setup_for_scoring_for_residue_opportunity( core::pose::Pose const & pose ) const;

void
--class--::setup_for_scoring_for_residue(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	ResSingleMinimizationData & min_data
) const;

bool
--class--::requires_a_setup_for_derivatives_for_residue_opportunity( core::pose::Pose const & pose ) const;

void
--class--::setup_for_derivatives_for_residue(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	ResSingleMinimizationData & min_data
) const;


bool
--class--::requires_a_setup_for_scoring_for_residue_pair_opportunity( core::pose::Pose const & pose ) const;

void
--class--::setup_for_scoring_for_residue_pair(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	ResSingleMinimizationData const & minsingle_data1,
	ResSingleMinimizationData const & minsingle_data2,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	ResPairMinimizationData & data_cache
) const;

bool
--class--::requires_a_setup_for_derivatives_for_residue_pair_opportunity( core::pose::Pose const & pose ) const;

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
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const & min_data,
	core::pose::Pose const & pose, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const;

void
--class--::backbone_backbone_energy(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const;

void
--class--::backbone_sidechain_energy(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const;

void
--class--::sidechain_sidechain_energy(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const;

bool
--class--::defines_intrares_energy( EnergyMap const & weights ) const = 0;

void
--class--::eval_intrares_energy(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const = 0;

bool
--class--::defines_intrares_energy_for_residue(
	core::conformation::Residue const & res
) const;

bool
--class--::use_extended_intrares_energy_interface() const;

void
--class--::eval_intrares_energy_ext(
	core::conformation::Residue const & rsd,
	ResSingleMinimizationData const & data_cache,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const;

void
--class--::eval_intrares_derivatives(
	core::conformation::Residue const & rsd,
	ResSingleMinimizationData const & min_data,
	core::pose::Pose const & pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const;

bool
--class--::defines_intrares_dof_derivatives( pose::Pose const & p ) const;

Real
--class--::eval_intraresidue_dof_derivative(
	core::conformation::Residue const & rsd,
	ResSingleMinimizationData const & min_data,
	core::id::DOF_ID const & dof_id,
	core::id::TorsionID const & torsion_id,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights
) const;

void
--class--::bump_energy_full(
	core::conformation::Residue const &,
	core::conformation::Residue const &,
	core::pose::Pose const &,
	ScoreFunction const &,
	EnergyMap &
) const;

virtual
void
--class--::bump_energy_backbone(
	core::conformation::Residue const &,
	core::conformation::Residue const &,
	core::pose::Pose const &,
	ScoreFunction const &,
	EnergyMap &
) const;

void
--class--::evaluate_rotamer_intrares_energies(
	core::conformation::RotamerSetBase const & set,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	utility::vector1< core::PackerEnergy > & energies
) const;

void
--class--::evaluate_rotamer_intrares_energy_maps(
	core::conformation::RotamerSetBase const & set,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	utility::vector1< EnergyMap > & emaps
) const;

void
--class--::evaluate_rotamer_pair_energies(
	core::conformation::RotamerSetBase const & set1,
	core::conformation::RotamerSetBase const & set2,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const;

void
--class--::evaluate_rotamer_background_energies(
	core::conformation::RotamerSetBase const & set,
	core::conformation::Residue const & residue,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const;

void
--class--::evaluate_rotamer_background_energy_maps(
	core::conformation::RotamerSetBase const & set,
	core::conformation::Residue const & residue,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	utility::vector1< EnergyMap > & emaps
) const;

--end_namespace--

