// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/util.hh
/// @brief Util functions for Input and Output.  Very general IO should go to utility/io.  These should be related to core in a deep way.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com), XRW 2016 Team

#ifndef INCLUDED_core_io_util_hh
#define INCLUDED_core_io_util_hh

// Project headers
#include <core/types.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/io/StructFileRep.fwd.hh>
#include <core/io/StructFileRepOptions.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Basic headers
#include <basic/Tracer.fwd.hh>

// Utility headers
#include <utility/io/ozstream.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <iosfwd>
#include <sstream>

namespace core {
namespace io {

///@brief Extract the pose energies table from an SFR as a string representation for PDB output.
std::string
pose_energies_from_sfr(
	StructFileRep const & sfr

);

///@brief Extract the pose energies table from an SFR as a string representation for PDB output.
void
pose_energies_from_sfr(
	StructFileRep const & sfr,
	std::stringstream & out
);


///@breif Extract the pose data cache from the SFR as a string representation for PDB output.
std::string
pose_data_cache_from_sfr(
	StructFileRep const & sfr
);

///@breif Extract the pose data cache from the SFR as a string representation for PDB output.
void
pose_data_cache_from_sfr(
	StructFileRep const & sfr,
	std::stringstream & out
);





void pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	utility::vector1< core::Size > const & residue_indices
);

void pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	utility::vector1< core::Size > const & residue_indices,
	StructFileRepOptions const & options
);

void pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< core::Size > const & residue_indices
);

void pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< core::Size > const & residue_indices,
	StructFileRepOptions const & options
);


} //core
} //io


#endif //core/io_util_hh

