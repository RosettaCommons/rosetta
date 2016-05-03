// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/PDBPoseOutputter.cc
/// @brief  Definition of the %PDBPoseOutputter class's methods
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

//unit headers
#include <protocols/jd3/pose_outputters/PDBPoseOutputter.hh>

//package headers
#include <protocols/jd3/LarvalJob.hh>

//project headers
#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>
//#include <core/io/pdb/build_pose_as_is.hh>

// ObjexxFCL
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <utility/io/ozstream.hh>

namespace protocols {
namespace jd3 {
namespace pose_outputters {

PDBPoseOutputter::PDBPoseOutputter() {}
PDBPoseOutputter::~PDBPoseOutputter() {}

bool PDBPoseOutputter::job_has_already_completed( LarvalJob const & /*job*/ ) const
{
	// STUBBED OUT!
	return false;
}


void PDBPoseOutputter::mark_job_as_having_started( LarvalJob const & /*job*/ ) const
{
	// STUBBED OUT!
}


void PDBPoseOutputter::write_output_pose( LarvalJob const & job, core::pose::Pose const & pose )
{
	std::string out_fname = output_pdb_name( job );
	utility::io::ozstream ostream( out_fname );

	/// Modified by VKM, 31 January 2016, Chemical XRW 2016.
	core::io::pdb::dump_pdb(
		pose,
		"",
		true,
		true,
		ostream,
		out_fname
	);

}

std::string
PDBPoseOutputter::output_pdb_name( LarvalJob const & job ) const
{
	return ( job.status_prefix() == "" ? "" : job.status_prefix() + "_" ) + job.job_tag() + "_" +
		ObjexxFCL::lead_zero_string_of( job.nstruct_index(), std::max( 4, int( std::log10( job.nstruct_max() ))) ) +
		".pdb";
}

} // namespace pose_outputters
} // namespace jd3
} // namespace protocols
