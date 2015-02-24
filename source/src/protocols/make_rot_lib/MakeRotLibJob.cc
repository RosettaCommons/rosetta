// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/make_rot_lib/MakeRotLibJob.fwd.hh
/// @brief  Implimentation file for MakeRotLibJob class
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

// unit headers
#include <protocols/make_rot_lib/MakeRotLibJob.hh>
#include <protocols/make_rot_lib/MakeRotLibOptionsData.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

namespace protocols {
namespace make_rot_lib {

MakeRotLibJob::MakeRotLibJob( jd2::InnerJobOP inner_job, core::Size nstruct_index,
	core::Real omg, utility::vector1< core::Real > bbs, utility::vector1< core::Size > bb_ids, core::Real eps,
	MakeRotLibOptionsDataOP mrlod, bool semirotameric ) :
	jd2::Job( inner_job, nstruct_index ),
	omg_( omg ),
	bbs_( bbs ),
	bb_ids_( bb_ids ),
	eps_( eps ),
	mrlod_( mrlod ),
	semirotameric_( semirotameric )
{}

}//make_rot_lib
}//protocols
