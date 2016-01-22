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
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperations.fwd.hh>

#include <basic/datacache/DataMap.hh>

#include <string>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>
#include <map>


namespace protocols {
namespace antibody {
namespace design {


/// @brief Get default options
AntibodyCDRSetOptions
get_cdr_set_options();

/// @brief Get options from an instruction file
AntibodyCDRSetOptions
get_cdr_set_options(std::string instruction_file);


/// @brief Get default options
AntibodyCDRGraftDesignOptions
get_graft_design_options();

/// @brief Get options from an instruction file
AntibodyCDRGraftDesignOptions
get_graft_design_options(std::string instruction_file);


/// @brief Get default options
AntibodyCDRSeqDesignOptions
get_seq_design_options();

/// @brief Get options from an instruction file
AntibodyCDRSeqDesignOptions
get_seq_design_options(std::string instruction_file);


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

/// @brief Disable design of the first 2 and last 3 residues of the H3 CDR (North CDR definitions - kink determinants)
core::pack::task::operation::RestrictResidueToRepackingOP
disable_h3_stem_positions(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose & pose,
	core::Size nter_stem = 2,
	core::Size cter_stem= 3);


core::scoring::ScoreFunctionOP
get_ab_design_global_scorefxn();

/// @brief Get Global Scorefxn from RS
core::scoring::ScoreFunctionOP
get_ab_design_global_scorefxn(utility::tag::TagCOP tag, basic::datacache::DataMap & data);

core::scoring::ScoreFunctionOP
get_ab_design_dock_high_scorefxn();

core::scoring::ScoreFunctionOP
get_ab_design_dock_low_scorefxn();

core::scoring::ScoreFunctionOP
get_ab_design_min_scorefxn();

core::scoring::ScoreFunctionOP
get_ab_design_min_scorefxn(utility::tag::TagCOP tag, basic::datacache::DataMap & data);

std::pair<bool, core::Size>
check_cb(core::pose::Pose const & pose, utility::vector1<bool> const & residues);

std::pair<bool, core::Size>
check_cb(core::pose::Pose const & pose, protocols::loops::Loops const & loops);

/// @brief Use insert_pose_into_pose to replace the cdr_piece with the current antibody's CDR.  No modeling or superposition.  For that, use protocols/grafting
void
insert_cdr_into_antibody(AntibodyInfoCOP ab_info, CDRNameEnum const cdr, core::pose::Pose & pose, core::pose::Pose & cdr_piece, core::Size overhang=3);


/// @brief Gets all possible graft permutations.
/// @details all_permutations is a list of vectors corresponding to cdrs_to_design vector.  Essentially, each inner index describes a position in the cdr_set.
/// Indexes correspond to CDRNameEnum, and whose values correspond to the cdr_set index.  If the value is 0, it means no cdr in set.
/// Example: <1, 0, 1, 1, 1, 1>.  This is a possible combination to try graft, the second CDR, H2 is not part of the combination.
utility::vector1<utility::vector1< core::Size > >
get_all_graft_permutations(

	utility::vector1<utility::vector1< core::Size > > permutations,
	utility::vector1<core::Size > totals,
	core::Size const n);

AntibodyDesignProtocolEnum
design_protocol_to_enum(std::string const & design_type);

std::string
design_protocol_to_string(AntibodyDesignProtocolEnum const design_type);

SeqDesignStrategyEnum
seq_design_strategy_to_enum(std::string const strategy);

std::string
seq_design_strategy_to_string(SeqDesignStrategyEnum strategy);

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

///These all need better names and they need to be moved to a general place.  Make the ResidueKey of PDBInfo a public class:

/// @brief Get the PDBNumbering from strings such as: 1A, 1A:I, 1:A:~, 1:A:I with A being chain A and I being an insertion code.
PDBNumbering
get_pdb_numbering_from_single_string( std::string const & pdb_residue);

/// @brief Get a resnum from strings such as: 1A, 1A:I, 1:A:~, 1:A:I with A being chain A and I being an insertion code.
core::Size
get_resnum_from_single_string(core::pose::Pose const & pose, std::string const & pdb_residue);

/// @brief Get a resnum using the PDBLandmark from strings such as: 1A, 1A:I, 1:A:~, 1:A:I with A being chain A and I being an insertion code.
core::Size
get_resnum_from_single_string_w_landmark(
	AntibodyInfoCOP ab_info,
	core::pose::Pose const & pose,
	std::string const & pdb_residue,
	AntibodyNumberingSchemeEnum const & scheme);

/// @brief Get a boolean vector of resnums with ranges, where a - indicates range.
/// Parses strings for PDB resnums such as 1A, 1A:I, 1:A:~, 1:A:I with A being chain A and I being an insertion code.
/// Example: 1A-10A or 1A-10A:I
///
utility::vector1<bool>
get_resnums_from_strings_with_ranges(core::pose::Pose const & pose, utility::vector1<std::string> const & pdb_residues);

/// @brief Get PDBNumbering from a vector of strings:
///  Example: 1A, 1A:I, 1:A:~, 1:A:I with A being chain A and I being an insertion code.
///
utility::vector1<PDBNumbering>
get_pdb_numbering_from_strings(utility::vector1<std::string> const & pdb_residues);

/// @brief get a boolean vector from a vector of strings:
///  Example: 1A, 1A:I, 1:A:~, 1:A:I with A being chain A and I being an insertion code.
///
utility::vector1<bool>
get_resnum_from_strings(core::pose::Pose const & pose, utility::vector1<std::string> const & pdb_residues);

/// @brief Get a boolean vector from a string of PDBNumbering.  No error checking.  Use with caution.
utility::vector1<bool>
get_resnum_from_pdb_numbering(core::pose::Pose const & pose, utility::vector1<PDBNumbering>const & pdb_residues);



///Add this across protocol...
//utility::vector1<PDBNumbering>
//get_pdb_numbering_from_strings_with_ranges(vector1<std::string> const & pdb_residues);

void
add_loops_from_bool_vector(loops::Loops & loops, utility::vector1< bool > residues, bool add_cutpoints = false);



/// @brief Get probability data for a given set of CDRs.  Will fill in the no_data_cdrs;
std::map< core::Size, std::map< core::chemical::AA, core::Real > >
get_cluster_profile_probability_data(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose& pose,
	utility::vector1<bool> const & cdrs,
	utility::vector1< bool > &no_data_cdrs,
	const core::Size prob_cutoff = 10,
	const bool use_outliers = false,
	const bool force_north_db = false);

//@brief Get probability data for a given set of CDRs.  Will fill in the no_data_cdrs;
std::map< core::Size, std::map< core::chemical::AA, core::Real > >
get_cluster_profile_probability_data(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose& pose,
	AntibodyCDRSeqDesignOptions const & seq_design_options,
	utility::vector1< bool > &no_data_cdrs,
	const core::Size prob_cutoff = 10,
	const bool use_outliers = false,
	const bool force_north_db = false);

/// @brief Transforms a sequence to a mutation set used by the AddCDRProfileSetsOperation.
/// Assumes that the sequence is the same length as the CDR.  Forces use of North CDR definitions.
std::map<core::Size, core::chemical::AA>
transform_sequence_to_mutation_set(
	AntibodyInfoCOP ab_info,
	core::pose::Pose const & pose,
	CDRNameEnum const cdr,
	std::string const & sequence);

/// @brief Set the native CDR sequence into the pose datacache
/// If none is set in the pose, will add it.
void
set_native_cdr_sequence(AntibodyInfoCOP ab_info, CDRNameEnum cdr, core::pose::Pose & pose);

std::string
get_native_sequence(core::pose::Pose const & pose);

bool
has_native_sequence(core::pose::Pose const & pose);

} //design
} //antibody
} //protocols


#endif //#ifndef INCLUDED_protocols/antibody/design_UTIL_HH
