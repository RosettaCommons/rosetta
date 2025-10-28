// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/util.hh
/// @brief Utility functions for the Antibody namespace
/// @author Jianqing Xu (xubest@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_util_hh
#define INCLUDED_protocols_antibody_util_hh


//Core Headers

#include <core/pose/Pose.fwd.hh>
#include <core/pose/DockingPartners.fwd.hh>
#include <core/types.hh>
#include <core/pack/task/TaskFactory.fwd.hh>


//Protocol Headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/AntibodyNumberingParser.fwd.hh>

//Utility Headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace antibody {



/// @brief Geta boolean vector of CDRs from a RosettaScripts tag. Companion function for XSD generation: attributes_for_get_cdr_bool_from_tag
utility::vector1<bool>
get_cdr_bool_from_tag(utility::tag::TagCOP tag, std::string const & name, bool include_cdr4 = false);

void
attributes_for_get_cdr_bool_from_tag(utility::tag::AttributeList& attlist,
	std::string const& tagname, std::string const& Description = "");

/// @brief Get a set of loops for a boolean vector of CDRNameEnums including any stem residues.
protocols::loops::LoopsOP
get_cdr_loops(
	AntibodyInfoCOP ab_info,
	core::pose::Pose const & pose,
	utility::vector1<bool> const & cdrs,
	core::Size stem_size = 0 );


////////////////////////////////////////// things to compare to native //////////////////////////////////////////
core::Real
global_loop_rmsd ( const core::pose::Pose & pose_in,
	const core::pose::Pose & native_pose,
	loops::LoopsOP current_loop );


/// @brief return false if any cdr cutpoint is broken
bool
cutpoints_separation( core::pose::Pose & pose, AntibodyInfoOP & antibody_info );

/// @details Compute the separation at the cutpoint. The N-C distance of the
///   peptide bond which should be formed at the cutpoint. A closed loop is
///   assumed to have a gap < 1.9 Ang
core::Real
cutpoint_separation(core::pose::Pose & pose_in, core::Size cutpoint);






/////////////////////////////////// Epitope + Paratope ////////////////////////////////////////////////////////////////////

/// @brief Get the epitope residues using the InterGroupNeighborsCalculator.
utility::vector1<bool>
select_epitope_residues(AntibodyInfoCOP ab_info, core::pose::Pose const & pose, core::Size const interface_distance = 10.0);




////////////////////////////////// Numbering ///////////////////////////////////

/// @brief Checks the length of the CDR, fixes AHO numbering for longer CDR loops that don't
/// fit within the numbering scheme
void
check_fix_aho_cdr_numbering(AntibodyInfoCOP ab_info, CDRNameEnum cdr, core::pose::Pose & pose);

void
check_fix_aho_cdr_numbering(AntibodyInfoCOP ab_info, core::pose::Pose & pose);







///////////////////////////////// Etc. /////////////////////////////////////
void
simple_one_loop_fold_tree(
	core::pose::Pose & pose,
	loops::Loop const & loop);

void
simple_fold_tree(
	core::pose::Pose & pose_in,
	core::Size jumppoint1,
	core::Size cutpoint,
	core::Size jumppoint2);

/// @brief Setup LH_A foldtree via docking.
core::pose::DockingPartners
setup_LH_A_foldtree(AntibodyInfoCOP ab_info, core::pose::Pose & pose);

/// @brief Setup A_LH foldtree via docking.
core::pose::DockingPartners
setup_A_LH_foldtree(AntibodyInfoCOP ab_info, core::pose::Pose & pose);


/// @brief Very specific packertask,
/// @details Ubound rotamer options, repack only, protein only, no disulfides.
core::pack::task::TaskFactoryOP
setup_packer_task( core::pose::Pose & pose_in);


/// @brief align current Fv to native.Fv
void
align_to_native( core::pose::Pose & pose,
	core::pose::Pose const & native_pose,
	AntibodyInfoOP const ab_info,
	AntibodyInfoOP const native_ab_info,
	std::string const & request_chain="LH");

bool
CDR_H3_filter_legacy_code_with_old_rule(
	const core::pose::Pose & pose_in,
	loops::Loop & input_loop,
	bool is_camelid);

bool
CDR_H3_cter_filter(
	const core::pose::Pose & pose_in,
	AntibodyInfoOP ab_info);

/// @brief Very basic kink function based on rama groups.  Kink is defined as having AB + DB at the C-terminal stem.
/// This is in the North/Dunbrack CDR definition.
/// Avoid this definition if you can and talk to brain.  This is currently used for a quick filtering for design.
bool
is_H3_rama_kinked(std::string const & rama);

/// @brief Apply auto-generated kink constraint
void
kink_constrain_antibody_H3( core::pose::Pose & pose, core::Size kink_begin );

void
kink_constrain_antibody_H3( core::pose::Pose & pose, AntibodyInfoOP const antibody_info );

void
qq_constrain_antibody( core::pose::Pose & pose, core::Size VH_qq_resi, core::Size VL_qq_resi );




///@brief get the equivalent landmark in a numbering scheme.  Set the resnum to zero if not found.
PDBLandmarkCOP
get_matching_landmark(
	AntibodyNumbering const & numbering,
	PDBLandmark const & landmark,
	AntibodyNumberingSchemeEnum const from_scheme,
	AntibodyNumberingSchemeEnum const to_scheme);

///@brief return vector of Chothia-nubmered, conserved residues for grafting and alignment
utility::vector1< core::Real >
get_conserved_residue_list( char chain );

} //namespace antibody
} //namespace protocols





#endif //INCLUDED_protocols_antibody_util_HH




