// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/util.hh
/// @brief  Utilities for JD2, which are somewhat safe to use outside of the JD2 system itself
/// @details For utilities internal to the jd2 system, see internal_util.hh
/// @author Andrew Leaver-Fay
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_jd2_util_hh
#define INCLUDED_protocols_jd2_util_hh

#include <core/pose/Pose.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/types.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <utility/vector1.hh>
#include <utility/file/FileName.hh>

#include <string>
#include <map>
#include <list>

namespace protocols {
namespace jd2 {

/// @brief is this application running with JD2
/// Useful for code that might not be running under JD2
bool jd2_used();

//////////////////////////////////////////////////////
//     JD2 transitional helpers
//
//  These are helper functions for accessing JD2 functionality in a potentially non-JD2 environment.
// They hopefully should be robust to being run in a non-JD2 environment.
//
//  Try not to use them if possible (but use them in preference to directly accessing JD2 objects)
//
//  If you do use them, use the fully qualified names (e.g. protocols::jd2::current_input_tag() ) so they are grep-able

/// @brief What is the tag of the current input structure?
std::string current_input_tag();

/// @brief What 'nstruct' value are we currently processing?
core::Size current_nstruct_index();

/// @brief What's the maximum nstruct we're processing?
core::Size max_nstruct_index();

/// @brief gets used output name of pose
/// ask jd for current-job ---> ask jobOutputter for name of this job
std::string current_output_name();

/// @brief call the 'filename' accessor of the current job-distributor with the current job
std::string current_output_filename();

std::string current_batch();

core::Size current_batch_id();

/// @brief returns 0 if no replicas (i.e., multiple processes per job )
/// otherwise it returns the sub-rank of the process within the job starting at 1
core::Size current_replica();

/// @brief Get the starting structure for the current run.
/// Will return a null pointer if that is not possible.
core::pose::PoseCOP get_current_jobs_starting_pose();

/// @brief writes pose to intermediate-scorefile using current Job and JobOutputter ..
/// @details copy count is used if multiple poses with same job-tag are written as for instance into a trajectory.
///  -1 indicates no copy count
///  >=0 copy_count will be attached as zerofilled postfix to job-tag
void output_intermediate_pose( core::pose::Pose const& pose, std::string const& stage_tag, int copy_count = -1, bool score_only = false );

/// @brief add output string
void add_string_to_current_job( std::string const & string_in );

/// @brief add output strings
void add_strings_to_current_job( std::list< std::string > const & strings );

/// @brief add a string/string pair
void add_string_string_pair_to_current_job( std::string const & string1, std::string const & string2 );

/// @brief add a string/real pair
void add_string_real_pair_to_current_job( std::string const & string_in, core::Real real_in );

std::list< std::string > get_strings_from_current_job();

std::map< std::string, std::string > get_string_string_pairs_from_current_job();

std::map< std::string, core::Real > get_string_real_pairs_from_current_job();

void
add_current_job_data_to_ss( core::io::silent::SilentStructOP ss );

void
write_score_tracer( core::pose::Pose const& pose_in, std::string tag );

///////////////////////////////////////
// Not actually JD2-specific

/// @brief Read the -s and -l flags to get a list of PDB files that should be
/// read in and operated upon.
utility::vector1< utility::file::FileName >
input_pdb_files_from_command_line();

void set_native_in_mover( protocols::moves::Mover &mover );

}  //jd2
}  //protocols

#endif
