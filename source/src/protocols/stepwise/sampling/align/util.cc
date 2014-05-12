// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/align/util.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/align/util.hh>
#include <protocols/stepwise/sampling/align/StepWisePoseAligner.hh>
#include <protocols/stepwise/sampling/util.hh>
#include <core/pose/Pose.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.sampling.align.util" );

namespace protocols {
namespace stepwise {
namespace sampling {
namespace align {

	///////////////////////////////////////////////////////////////////////
	// By using the StepWisePoseAligner,
	// this is really careful with changes in virtual atoms, etc.
	///////////////////////////////////////////////////////////////////////
	core::Real
	get_rmsd( core::pose::Pose const & pose1, core::pose::Pose const & pose2,
						utility::vector1< core::Size > const & calc_rms_res,
						bool const check_align_at_superimpose_res,
						bool const check_switch ) {

		StepWisePoseAligner pose_aligner( pose2 );
		pose_aligner.set_user_defined_calc_rms_res( calc_rms_res );
		if ( check_align_at_superimpose_res )	{
			pose_aligner.set_root_partition_res( sampling::figure_out_root_partition_res( pose2, calc_rms_res ) );
		}
		pose_aligner.initialize( pose1 );
		core::Real const rmsd = pose_aligner.get_rmsd_no_superimpose( pose1, check_align_at_superimpose_res );

		if ( check_switch ) {
			core::Real const rmsd_switch = get_rmsd( pose2, pose1, calc_rms_res,
																							 check_align_at_superimpose_res, false /* check switch */);
			runtime_assert( std::abs( rmsd  - rmsd_switch ) < 1.0e-3 );
		}

		return rmsd;
	}

	///////////////////////////////////////////////////////////////////////
	core::Real
	get_rmsd( core::pose::Pose const & pose1, core::pose::Pose const & pose2,
						bool const check_align_at_superimpose_res,
						bool const check_switch ) {
		utility::vector1< core::Size > blank_calc_rms_res;
		return get_rmsd( pose1, pose2, blank_calc_rms_res,
										 check_align_at_superimpose_res, check_switch );
	}

} //align
} //sampling
} //stepwise
} //protocols
