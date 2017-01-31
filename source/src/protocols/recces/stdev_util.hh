// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/stdev_util.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_recces_stdev_util_HH
#define INCLUDED_protocols_recces_stdev_util_HH

#include <protocols/recces/sampler/MC_Comb.fwd.hh>
#include <protocols/recces/options/RECCES_Options.fwd.hh>
#include <protocols/recces/params/RECCES_Parameters.fwd.hh>
#include <protocols/recces/sampler/MC_Sampler.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/rna/RNA_SecStruct.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace recces {

void
set_sampler_gaussian_stdev( protocols::recces::sampler::MC_CombOP sampler,
	core::Real const & temperature,
	core::pose::Pose const & pose,
	protocols::recces::options::RECCES_Options const & options,
	protocols::recces::params::RECCES_Parameters const & params );


void
set_gaussian_stdevs_recces_turner( protocols::recces::sampler::MC_CombOP sampler,
	core::Real const & temperature,
	core::pose::Pose const & pose,
	core::pose::rna::RNA_SecStruct const & rna_secstruct,
	protocols::recces::params::RECCES_Parameters const & params );

core::Real gaussian_stdev( core::Real const n_rsd, core::Real const temp, bool const is_bp );

void
set_gaussian_stdevs_thermal_sampler(
	protocols::recces::sampler::MC_SamplerOP sampler,
	core::Real const & temperature,
	core::pose::Pose const & pose,
	protocols::recces::options::RECCES_Options const & options );

} //recces
} //protocols

#endif
