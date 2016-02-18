// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/carbohydrates/pose_io.hh
/// @brief   Pose input/output function declarations for carbohydrate-specific data formats.
/// @author  Labonte <JWLabonte@jhu.edu>

#ifndef INCLUDED_core_io_carbohydrates_pose_io_HH
#define INCLUDED_core_io_carbohydrates_pose_io_HH

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.hh>


namespace core {
namespace io {
namespace carbohydrates {

// Input //////////////////////////////////////////////////////////////////////

// TODO: Add file input methods below.


/// @brief Parse sugar code suffixes to extract a list of sugar modifications with their corresponding positions.
utility::vector1< std::pair< core::uint, std::string > > sugar_modifications_from_suffix( std::string const & suffix );


/// @brief  Read a CFG structure from a GWS-formatted string and load into a pose.
//void pose_from_gws_string( core::pose::Pose & pose, std::string const & filename );

/// @brief  Read the first CFG structure from a GWS file and load into a pose.
//void pose_from_gws( core::pose::Pose & pose, std::string const & filename );

/// @brief  Read all CFG structures from a GWS file and load into poses.
//void poses_from_gws( utility::vector1< core::pose::Pose > & poses, std::string const & filename );

/// @brief  Create a pose from a GWS file.
//core::pose::PoseOP pose_from_gws( std::string const & filename );


// Output /////////////////////////////////////////////////////////////////////

/// @brief  Return a GWS-formatted string for the given carbohydrate residue.
std::string residue_gws_string( core::pose::Pose const & pose, core::uint const seqpos );

/// @brief  Return a GWS-formatted string for each carbohydrate residue in the given sequence range, including branches.
std::string residue_range_gws_string( core::pose::Pose const & pose, core::uint const begin, core::uint const end );

/// @brief  Return a GWS-formatted string for the given carbohydrate chain, including branches.
std::string chain_gws_string( core::pose::Pose const & pose, core::uint const chain_id );


/// @brief  Write the GlycoWorkbench structure file for the given pose chain to <filename>.
void dump_gws_chain( core::pose::Pose const & pose, core::uint const chain_id, std::string const & filename );

/// @brief  Write the GlycoWorkbench structure file for all carbohydrate chains of the given pose to <filename>.
void dump_gws( core::pose::Pose const & pose, std::string const & filename );

// Utility /////////////////////////////////////////////////////////////////////

/// @brief Given a char, parse it as an integer.
/// @details Returns 0 for anything outside of the range 1-9.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
core::uint char_to_int( char const char_in );

}  // namespace carbohydrates
}  // namespace io
}  // namespace core

#endif  // INCLUDED_core_io_carbohydrates_pose_io_HH
