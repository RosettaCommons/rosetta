// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/internal_util.hh
/// @brief  Utilities for JD2, intended to be used within the JD2 system itself
/// @details For utilities which might be used outside JD2, see util.hh

#ifndef INCLUDED_protocols_jd2_internal_util_hh
#define INCLUDED_protocols_jd2_internal_util_hh

#ifdef USEMPI
#include <mpi.h>
#endif

#include <core/pose/Pose.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/types.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <utility/vector1.hh>
#include <utility/file/FileName.hh>

#include <string>
#include <map>
#include <list>

namespace protocols {
namespace jd2 {

//////////////////////////////////////////////////////
//
//  Please avoid using these outside of JD2-specific code


void add_job_data_to_ss( core::io::silent::SilentStructOP ss, JobCOP job_op );

void register_options();

JobOP get_current_job();

#ifdef USEMPI
/// @brief returns communicator defined by the JobDistributor or MPI_COMM_WORLD
MPI_Comm const& current_mpi_comm();
#endif

}  //jd2
}  //protocols

#endif
