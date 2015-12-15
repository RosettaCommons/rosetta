// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/pose_inputters/PDBPoseInputter.hh
/// @brief  Declaration of the %PDBPoseInputter class for initializing Poses from .pdb or .pdb.gz files
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Steven Lewis (smlewi@gmail.com)


#ifndef INCLUDED_protocols_jd3_pose_inputters_PDBPoseInputter_HH
#define INCLUDED_protocols_jd3_pose_inputters_PDBPoseInputter_HH

// Unit headers
#include <protocols/jd3/pose_inputters/PDBPoseInputter.fwd.hh>

// Package headers
#include <protocols/jd3/PoseInputter.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

//utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace jd3 {
namespace pose_inputters {

/// @brief This is the simplest implementation of PoseInputter, which reads from -s/-l and PDB files.
class PDBPoseInputter : public protocols::jd3::PoseInputter
{
public:

	PDBPoseInputter();
	virtual ~PDBPoseInputter();

	/// @brief Constructs a list of PoseInputSource objects reading from the
	/// -s or -l command line flags. This stores the names of the PDBs that
	/// are to be read in, and it initializes the input tags based on the pdb
	/// names, stripping the path and the extension from the file name.
	PoseInputSources initialize_pose_input_sources();

	/// @brief Takes a PoseInputSource object previously initialized in the
	/// call to initialize_pose_input_sources()
	core::pose::PoseOP pose_from_input_source( PoseInputSource const & );

}; // PDBPoseInputter

} // namespace pose_inputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_PDBPoseInputter_HH
