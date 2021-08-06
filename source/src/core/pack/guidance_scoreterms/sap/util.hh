// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/guidance_scoreterms/sap/util.hh
/// @brief  Utility functions for SapScore
/// @author Brian Coventry (bcov@uw.edu)


#ifndef INCLUDED_core_pack_guidance_scoreterms_sap_util_hh
#define INCLUDED_core_pack_guidance_scoreterms_sap_util_hh

// Unit headers
#include <core/pack/guidance_scoreterms/sap/SapConstraintOptions.fwd.hh>
#include <core/pack/guidance_scoreterms/sap/SapParameterOptions.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintHelper.fwd.hh>

// Package headers
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

#include <core/conformation/Residue.fwd.hh> // AUTO IWYU For ResidueCOP
#include <core/scoring/annealing/RotamerSets.fwd.hh> // AUTO IWYU For RotamerSetsOP



namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {

pack::rotamer_set::RotamerSetsOP
rotamer_sets_from_pose(
	pose::Pose const & pose,
	utility::vector1< core::conformation::ResidueCOP > & res_vector
);

SapConstraintHelperOP
common_setup(
	SapConstraintOptionsOP const & options,
	pose::Pose const & pose,
	select::residue_selector::ResidueSelectorCOP const & score_sel,
	select::residue_selector::ResidueSelectorCOP const & sap_calculate_sel,
	select::residue_selector::ResidueSelectorCOP const & sasa_sel,
	utility::vector1< core::conformation::ResidueCOP > & res_vector
);


Real
calculate_sap(
	pose::Pose const & pose,
	select::residue_selector::ResidueSelectorCOP const & score_sel,
	select::residue_selector::ResidueSelectorCOP const & sap_calculate_sel,
	select::residue_selector::ResidueSelectorCOP const & sasa_sel,
	SapParameterOptions const & sap_parameter_options = SapParameterOptions()
);

core::id::AtomID_Map<Real>
calculate_per_atom_sap(
	pose::Pose const & pose,
	select::residue_selector::ResidueSelectorCOP const & score_sel,
	select::residue_selector::ResidueSelectorCOP const & sap_calculate_sel,
	select::residue_selector::ResidueSelectorCOP const & sasa_sel,
	SapParameterOptions const & sap_parameter_options = SapParameterOptions()
);

utility::vector1<Real>
calculate_per_res_sap(
	pose::Pose const & pose,
	select::residue_selector::ResidueSelectorCOP const & score_sel,
	select::residue_selector::ResidueSelectorCOP const & sap_calculate_sel,
	select::residue_selector::ResidueSelectorCOP const & sasa_sel,
	SapParameterOptions const & sap_parameter_options = SapParameterOptions()
);

Real
calculate_slow_approx_sap(
	pose::Pose const & pose,
	select::residue_selector::ResidueSelectorCOP const & score_sel,
	select::residue_selector::ResidueSelectorCOP const & sap_calculate_sel,
	select::residue_selector::ResidueSelectorCOP const & sasa_sel,
	SapParameterOptions const & sap_parameter_options = SapParameterOptions()
);

core::id::AtomID_Map<Real>
sap_atom_sasa(
	core::pose::Pose const & pose,
	select::residue_selector::ResidueSubset const & sasa_sub
);

} //sap
} //guidance_scoreterms
} //pack
} //core

#endif
