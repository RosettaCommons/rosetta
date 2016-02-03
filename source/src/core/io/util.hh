// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

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

///@brief specific precision for score output.  Don't use this, use string_util, fmt_real or Real2string
core::Real restrict_prec( core::Real inval );

/// @brief Write extra Pose energies information into a string and return it.
/// @details Added during the 2016 Chemical XRW.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
std::string extract_extra_scores(
	pose::Pose const & pose
);

void extract_extra_scores(
	core::pose::Pose const & pose,
	std::stringstream & out
);


/// @brief Write Pose energies information into a string and return it.
/// @details Added during the 2016 Chemical XRW.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
std::string extract_scores(
	core::pose::Pose const & pose,
	std::string const &filename=""
);

/// @brief Write  <pose>  Energies information into an output stream
/// (e.g. the tail of a pdb file)
void extract_scores(
	core::pose::Pose const & pose,
	std::stringstream & out,
	std::string const &filename=""
);

/// @brief Utility function to round a real value to the given precisions (number of digits after the decimal place) for output.
/// For use solely by extract_scores()
/// @details Apparently, there isn't an easy way to do this in C++, or even the general goal
/// of limiting the precision in output streams. (setprecision() with default formatting doesn't
/// correctly handle very small numbers, and with fixed precision outputs superfluous zeros.)
///
/// Do Not use this!  See utility::string_util, fmt_real and Real2string for precision output.
core::Real
restrict_prec( core::Real inval );

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

