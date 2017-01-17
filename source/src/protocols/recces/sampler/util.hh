// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/sampler_util.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_recces_sampler_util_HH
#define INCLUDED_protocols_recces_sampler_util_HH

#include <protocols/recces/params/RECCES_Parameters.fwd.hh>
#include <protocols/recces/options/RECCES_Options.fwd.hh>
#include <protocols/recces/sampler/MC_Comb.fwd.hh>
#include <protocols/recces/sampler/rna/MC_RNA_OneJump.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/rna/RNA_SecStruct.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace recces {
namespace sampler {

///////////////////////////////////////////////////////////////////////////////////////////////////
// Sampler setup
protocols::recces::sampler::MC_CombOP
initialize_sampler( core::pose::Pose const & pose,
										options::RECCES_Options const & options,
										params::RECCES_Parameters const & params );

protocols::recces::sampler::MC_CombOP
get_recces_turner_sampler( core::pose::Pose const & pose,
													 core::Real const & a_form_range,
													 core::pose::rna::RNA_SecStruct const & secstruct,
													 params::RECCES_Parameters const & params );

protocols::recces::sampler::MC_CombOP
initialize_thermal_sampler( core::pose::Pose const & pose,
	options::RECCES_Options const & options );

protocols::recces::sampler::rna::MC_RNA_OneJumpOP
initialize_jump_sampler( core::pose::Pose const & pose,
												 core::Size const & num_jump,
												 options::RECCES_Options const & options );

} //sampler
} //recces
} //protocols

#endif
