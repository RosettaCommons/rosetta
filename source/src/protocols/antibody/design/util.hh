// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/design/util.hh
/// @brief Utility functions for antibody design namespace.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_design_UTIL_HH
#define INCLUDED_protocols_antibody_design_UTIL_HH

#include <protocols/antibody/design/AntibodyDesignEnum.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/database/CDRSetOptions.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.hh>
#include <protocols/antibody/design/CDRGraftDesignOptions.hh>
#include <protocols/loops/Loops.hh>

#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperations.fwd.hh>

#include <string>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <map>


namespace protocols {
namespace antibody{
namespace design{


/// @brief Use insert_pose_into_pose to replace the cdr_piece with the current antibody's CDR.  No modeling or superposition.  For that, use protocols/grafting
void
insert_cdr_into_antibody(AntibodyInfoCOP ab_info, CDRNameEnum const cdr, core::pose::Pose & pose, core::pose::Pose & cdr_piece, core::Size overhang=3);


/// @brief Gets all possible graft permutations.
/// @details all_permutations is a list of vectors corresponding to cdrs_to_design vector.  Essentially, each inner index describes a position in the cdr_set.
/// Indexes correspond to CDRNameEnum, and whose values correspond to the cdr_set index.  If the value is 0, it means no cdr in set.
/// Example: <1, 0, 1, 1, 1, 1>.  This is a possible combination to try graft, the second CDR, H2 is not part of the combination.
void
get_all_graft_permutations(
	utility::vector1<core::Size > & total_cdr_set,
	utility::vector1<vector1< core::Size > > & all_permutations,
	utility::vector1< core::Size >current_index,
	core::Size const cdr_num);

DesignTypeEnum
design_type_from_string(std::string const design_type);

// Undefined, commenting out to fix PyRosetta build  std::string design_type_from_enum(DesignTypeEnum const design_type);

/// @brief Convert an ab_dock_chain (L_H/ LH_A, etc. to the full dock chain string)
std::string
get_dock_chains_from_ab_dock_chains(AntibodyInfoCOP ab_info, std::string ab_dock_chains);


/// Move this somewhere in pose or pose_selection.  PDBInfo's ResidueKey should be public and passed around in the first place.
/// @brief Super-basic numbering struct.
struct PDBNumbering {
	core::Size resnum;
	char chain;
	char icode;
};

//Get PDBNumbering from a string that looks like 10A for resnum and chain; 10A:A for optional icode.
utility::vector1<PDBNumbering>
get_pdb_numbering_from_string(vector1<std::string> const & pdb_residues);

///Get a boolean vector from a string of PDBNumbering.  No error checking.  Use with caution.
utility::vector1<bool>
get_resnum_from_pdb_numbering(core::pose::Pose const & pose, vector1<PDBNumbering>const & pdb_residues);

std::pair<bool, core::Size>
check_cb(core::pose::Pose const & pose, utility::vector1<bool> const & residues);

std::pair<bool, core::Size>
check_cb(core::pose::Pose const & pose, protocols::loops::Loops const & loops);

/// @brief Disable design of any particular region of the antibody complex.
core::pack::task::operation::RestrictResidueToRepackingOP
disable_design_region(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose & pose,
	AntibodyRegionEnum region);

/// @brief Disable design of the antigen residues
core::pack::task::operation::RestrictResidueToRepackingOP
disable_design_antigen(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose & pose);

/// @brief Disable design of the framework residues
core::pack::task::operation::RestrictResidueToRepackingOP
disable_design_framework(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose & pose);

/// @brief Get a Restrict operation to turn OFF design for all CDRs.
core::pack::task::operation::RestrictResidueToRepackingOP
disable_design_cdrs(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose & pose);

/// @brief Get a Restrict operation to turn OFF design for particular CDRs.
core::pack::task::operation::RestrictResidueToRepackingOP
disable_design_cdr(
	AntibodyInfoCOP ab_info,
	CDRNameEnum cdr,
	const core::pose::Pose & pose);

/// @brief Disable design for conserved framework positions.
/// TODO: Needs to be expanded to read positions from database.
core::pack::task::operation::RestrictResidueToRepackingOP
disable_conserved_framework_positions(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose & pose);

/// Application Options - Should be moved as part of parsers? 

/// @brief Get options set from default instructions file and any user overrides
AntibodyCDRSetOptions
get_cdr_set_options();

/// @brief Get options set from default instructions file and any user overrides
AntibodyCDRGraftDesignOptions
get_graft_design_options();

/// @brief Get options set from default instructions file and any user overrides
AntibodyCDRSeqDesignOptions
get_seq_design_options();


} //design
} //antibody
} //protocols


#endif	//#ifndef INCLUDED_protocols/antibody/design_UTIL_HH
