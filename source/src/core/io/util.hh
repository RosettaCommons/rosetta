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

/// @brief Identify residues that are branch points (have more than one downstream neighbor)
/// One eneighbpr will be considered as the mainchain continuations. All others are branches.
void
find_branch_points( Size const & seqpos, chemical::ResidueTypeCOP & RT, bool & is_branch_point, utility::vector1< std::string > & branch_points_on_this_residue, utility::vector1< std::string > const & rosetta_residue_name3s, Size mainchain_neighbor, utility::vector1< core::io::ResidueInformation > & rinfos_, StructFileRep::Strings& branch_lower_termini_, utility::vector1< Size >& glycan_positions_, core::io::StructFileRepOptions const & options, StructFileRep const & sfr );

/// @brief Helper function to find connected residues.
void
find_downstream_neighbor( core::Size const seqpos, Vector const & upstream_atom_xyz, std::pair<core::Size, std::string> &neighbor, utility::vector1< core::io::ResidueInformation > const & rinfos_, utility::vector1< Size > const & glycan_positions_, utility::vector1< std::string > const & rosetta_residue_name3s, core::io::StructFileRepOptions const & options );

/// @detail Loop through a residue's atoms to see if there are othe residues connected to it.
/// Then the residue position that contiues the main chain (if any) is returned.
Size
find_mainchain_connection( utility::vector1< core::io::ResidueInformation >& rinfos_, core::io::StructFileRep& sfr_,  std::string const & resid, Size const & seqpos, utility::vector1< std::string > const & rosetta_residue_name3s, bool & same_chain_next, bool & is_upper_terminus, int const CARB_MAINCHAIN_CONN_POS, utility::vector1< Size >& glycan_positions_, core::io::StructFileRepOptions const & options );

/// @brief Recursively find a child residue and it's children and it's children ....
/// This function figures out in which order glycan residues are connected
void
find_children( Size const seqpos, utility::vector1< core::io::ResidueInformation > const & rinfos_, chemical::ResidueTypeSetCOP residue_type_set_, utility::vector1< std::string > const & rosetta_residue_name3s_, utility::vector1< core::Size > & correct_order_, utility::vector1< Size > const & glycan_positions_, utility::vector1< core::Size > & glycan_positions_temp_, core::io::StructFileRepOptions const & options );

/// @brief Test if a glycan residue is the first of it's tree
/// This is done by checking if there is anything conncted at the C1 position
bool is_root( core::Size const seqpos, utility::vector1< core::io::ResidueInformation > const & rinfos_, utility::vector1< Size > const & glycan_positions_, StructFileRep::Strings & branch_lower_termini_extra_, core::io::StructFileRepOptions const & options );

/// @brief Function to identify glycans and fix their rosetta names
/// This function determines the correct name3s and fills in the glycan_positions_ vector
void
fix_residue_info_and_order( utility::vector1< core::io::ResidueInformation >& rinfos, core::io::StructFileRep& sfr, chemical::ResidueTypeSetCOP residue_type_set, utility::vector1< std::string >& rosetta_residue_name3s, core::io::StructFileRep::Strings & branch_lower_termini_extra, utility::vector1< std::string > & glycan_tree_roots, utility::vector1< core::Size >& glycan_positions, core::io::StructFileRepOptions const & options );
/// @brief Bring glycans into the correct order, which corresponds to connectivity of ech glycan tree
/// This requires reordering rinfos_ and rosetta_residue_name3s_.
void reorder_glycan_residues( utility::vector1< core::io::ResidueInformation >& rinfos_, utility::vector1< std::string >& rosetta_residue_name3s_, utility::vector1< core::Size >& correct_order_, utility::vector1< core::Size > const & glycan_positions_  );

/// @brief sort link records in srf.link_map()
void sort_link_records( utility::vector1<core::io::LinkInformation> & link_records );

} //core
} //io


#endif //core/io_util_hh

