// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/util.hh
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Base class
/// @author Andrew Leaver-Fay
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_jd2_util_hh
#define INCLUDED_protocols_jd2_util_hh

#ifdef USEMPI
#include <mpi.h>
#endif

#include <core/pose/Pose.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/types.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <utility/vector1.hh>


#ifdef WIN32
#include <string>
#endif


namespace protocols {
namespace jd2 {

/// @brief writes pose to intermediate-scorefile using current Job and JobOutputter ..
/// @detail copy count is used if multiple poses with same job-tag are written as for instance into a trajectory.
///  -1 indicates no copy count
///  >=0 copy_count will be attached as zerofilled postfix to job-tag
void output_intermediate_pose( core::pose::Pose const& pose, std::string const& stage_tag, int copy_count = -1, bool score_only = false );

/// @brief gets used output name of pose
/// ask jd for current-job ---> ask jobOutputter for name of this job
std::string current_output_name();

/// @brief call the 'filename' accessor of the current job-distributor with the current job
std::string current_output_filename();

std::string current_batch();

/// @brief is this application running with jd2 --- used for some code that yields backward compatability with old JobDistributor
bool jd2_used();

core::pose::PoseCOP get_current_jobs_starting_pose();

void add_job_data_to_ss( core::io::silent::SilentStructOP ss, JobCOP job_op );

void register_options();

JobOP get_current_job();


void set_native_in_mover( protocols::moves::Mover &mover );

void
write_score_tracer( core::pose::Pose const& pose_in, std::string tag );

#ifdef USEMPI
/// @brief returns communicator defined by the JobDistributor or MPI_COMM_WORLD
MPI_Comm const& current_mpi_comm();
#endif

/// @brief returns 0 if no replicas (i.e., multiple processes per job )
/// otherwise it returns the sub-rank of the process within the job starting at 1
core::Size current_replica();

}  //jd2
}  //protocols

#endif
