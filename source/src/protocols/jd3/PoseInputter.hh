// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/PoseInputter.hh
/// @brief  Declaration of the %PoseInputter class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_PoseInputter_HH
#define INCLUDED_protocols_jd3_PoseInputter_HH

//unit headers
#include <protocols/jd3/PoseInputter.fwd.hh>

// Package headers
#include <protocols/jd3/PoseInputSource.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

// C++ headers
#include <map>

namespace protocols {
namespace jd3 {

/// @brief The %PoseInputter is responsible for reading from the command line a set of structures
/// that are each to be run through a protocol, where each input struture will be the starting
/// point for some number of jobs (where that number is at the JobQueen's discretion).  The
/// %PoseInputter is responsible for two things:
/// - for creating a list of PoseInputSource objects
/// - for turning a PoseInputSource object into a Pose on demand.
class PoseInputter
{
public:

	PoseInputter();
	virtual ~PoseInputter();

	/// @brief Construct a list of PoseInputSource objects that will be used by the
	/// JobQueen to, from there, construct a list of Jobs.  This will invariably
	/// rely upon the command line.
	virtual PoseInputSources initialize_pose_input_sources() = 0;

	/// @brief Convert a single PoseInputSource into a Pose that will be used to
	/// initialize a Job.  The PoseInputSource object must have originated from
	/// this %PoseInputter in the prior call to initialize_pose_input_sources.
	virtual core::pose::PoseOP pose_from_input_source( PoseInputSource const & ) = 0;

}; // PoseInputter

} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_PoseInputter_HH
