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
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRep.fwd.hh>
#include <core/io/StructFileRepOptions.fwd.hh>
#include <core/io/ResidueInformation.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
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

/////////////////////////////////////////

/// @brief Reorganize the given glycan residues such that they're in tree order
/// @details Not 100% great, as it won't move any non-glycan residue, which may result in
/// a funky order if there are non-glycans interspersed with the glycan residues
/// Return the positions (in the rinfos) of the known (non-reducing) ends of the glycan chains
utility::vector1< core::Size >
fix_glycan_order( utility::vector1< core::io::ResidueInformation > & rinfos,
	utility::vector1< core::Size > const & glycan_positions,
	StructFileRepOptions const & options,
	std::map< std::string, std::map< std::string, std::pair< std::string, std::string > > > known_links );

/// @brief Find the order of glycan residues for all glycan positions
utility::vector1< core::Size >
find_carbohydrate_order(
	utility::vector1< core::io::ResidueInformation > const & rinfos,
	utility::vector1< core::Size > const & glycan_positions,
	utility::vector1< core::Size > & chain_ends, // return-by-reference of (rinfos) positions of (non-reducing) end sugars
	// map of anomeric positions to where they are connected to
	std::map< std::pair< core::Size, std::string >, std::pair< core::Size, std::string > > const & link_map,
	std::map< std::string, std::map< std::string, std::pair< std::string, std::string > > > known_links );

/// @brief Find the order of all the glycans connected to current_res, as indicated by the connectivity map.
utility::vector1< core::Size >
find_carbohydrate_subbranch_order( core::Size current_res,
	utility::vector1< core::Size > & chain_ends, // return-by-reference for (non-reducing) end sugars
	// Nested map positions & (non-anomeric) positions to attached anomeric position
	std::map< core::Size, std::map< std::string, std::pair< core::Size, std::string > > > const & connectivity,
	std::set< core::Size > & addressed );

/// @brief Bring glycans into the correct order, which corresponds to connectivity of each glycan tree
/// This requires reordering rinfos and rosetta_residue_name3s.
void reorder_glycan_residues( utility::vector1< core::io::ResidueInformation > & rinfos,
	utility::vector1< core::Size > const & correct_order,
	utility::vector1< core::Size > const & glycan_positions );

/// @brief Determine links between glycan residues based on coordinates
/// Returns a map keyed on anomeric positions to the atom it's nominally attached to
std::map< std::pair< core::Size, std::string >, std::pair< core::Size, std::string > >
determine_glycan_links( utility::vector1< core::io::ResidueInformation > const & rinfos,
	StructFileRepOptions const & options );

std::map< std::string, std::map< std::string, std::pair< std::string, std::string > > >
explicit_links_from_sfr_linkage( std::map< std::string, utility::vector1< LinkInformation > > const & link_map,
	utility::vector1< core::io::ResidueInformation > const & rinfos );

void
add_glycan_links_to_map(
	std::map< std::string, std::map< std::string, std::pair< std::string, std::string > > > & known_links,
	std::map< std::pair< core::Size, std::string >, std::pair< core::Size, std::string > > const & link_map,
	utility::vector1< core::io::ResidueInformation > const & rinfos );

} //core
} //io


#endif //core/io_util_hh

