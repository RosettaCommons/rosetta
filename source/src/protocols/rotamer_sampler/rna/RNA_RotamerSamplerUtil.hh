// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_RotamerSamplerUtil.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rotamer_sampler_rna_RNA_RotamerSamplerUtil_HH
#define INCLUDED_protocols_rotamer_sampler_rna_RNA_RotamerSamplerUtil_HH

#include <protocols/rotamer_sampler/rna/RNA_RotamerSamplerUtil.fwd.hh>
#include <protocols/rotamer_sampler/RotamerBase.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_ModelerOptions.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace rotamer_sampler {
namespace rna {

	rotamer_sampler::RotamerBaseOP
	setup_rotamer_sampler( core::pose::Pose const & pose,
												 protocols::stepwise::sampling::rna::StepWiseRNA_ModelerOptionsCOP options,
												 protocols::stepwise::sampling::rna::StepWiseRNA_JobParametersCOP job_parameters,
												 bool const build_pose_from_scratch,
												 bool const kic_sampling,
												 bool const close_chain );

} //rna
} //rotamer_sampler
} //protocols

#endif
