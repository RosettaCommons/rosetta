// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

#ifndef INCLUDED_protocols_recces_thermal_sampler_HH
#define INCLUDED_protocols_recces_thermal_sampler_HH

// libRosetta headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>

#include <protocols/recces/sampler/rna/MC_RNA_MultiSuite.hh>
#include <protocols/moves/SimulatedTempering.hh>
#include <protocols/recces/sampler/rna/MC_RNA_KIC_Sampler.hh>

#include <protocols/recces/util.hh>
#include <protocols/recces/sampler/MC_OneTorsion.hh>

// C++ headers
#include <iostream>
#include <string>

namespace protocols {
namespace recces {

//////////////////////////////////////////////////////////////////////////////
utility::vector1<core::Real> get_torsions(
	utility::vector1<core::id::TorsionID> & torsion_ids,
	const core::pose::Pose & pose
);

//////////////////////////////////////////////////////////////////////////////
void set_gaussian_stdevs(
	utility::vector1<protocols::recces::sampler::rna::MC_RNA_KIC_SamplerOP> & internal_bb_sampler,
	utility::vector1<protocols::recces::sampler::MC_OneTorsionOP> & chi_sampler,
	protocols::recces::sampler::rna::MC_RNA_MultiSuite & standard_bb_sampler,
	moves::SimulatedTempering const & tempering,
	core::Size const & total_rsd,
	core::Size const & sampled_rsd,
	utility::vector1<bool> is_free
);

} //recces
} //protocols

#endif

