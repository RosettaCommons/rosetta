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

#ifndef INCLUDED_protocols_farna_thermal_sampling_thermal_sampler_HH
#define INCLUDED_protocols_farna_thermal_sampling_thermal_sampler_HH

// libRosetta headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>

#include <protocols/stepwise/sampler/rna/RNA_MC_MultiSuite.hh>
#include <protocols/moves/SimulatedTempering.hh>
#include <protocols/stepwise/sampler/rna/RNA_MC_KIC_Sampler.hh>

#include <protocols/farna/thermal_sampling/util.hh>
#include <protocols/stepwise/sampler/MC_OneTorsion.hh>

// C++ headers
#include <iostream>
#include <string>

namespace protocols {
namespace farna {
namespace thermal_sampling {

using namespace core::pose;
using namespace core;
using namespace protocols;
using namespace protocols::stepwise;

//////////////////////////////////////////////////////////////////////////////
utility::vector1<core::Real> get_torsions(
	utility::vector1<core::id::TorsionID> & torsion_ids,
	const Pose & pose
);

//////////////////////////////////////////////////////////////////////////////
void set_gaussian_stdevs(
	utility::vector1<protocols::stepwise::sampler::rna::RNA_MC_KIC_SamplerOP> & internal_bb_sampler,
	utility::vector1<protocols::stepwise::sampler::MC_OneTorsionOP> & chi_sampler,
	sampler::rna::RNA_MC_MultiSuite & standard_bb_sampler,
	moves::SimulatedTempering const & tempering,
	Size const & total_rsd,
	Size const & sampled_rsd,
	utility::vector1<bool> is_free
);

}
}
}

#endif

