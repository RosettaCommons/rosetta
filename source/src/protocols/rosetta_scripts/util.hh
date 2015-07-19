// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/RosettaScripts/util.hh
/// @brief Utility functions useful in RosettaScripts.
/// @author Sarel Fleishman (sarelf@u.washington.edu)
///         Jacob Corn (jecorn@u.washington.edu)
///         Rocco Moretti (rmoretti@u.washington.edu)
///         Eva-Maria Strauch (evas01@uw.edu)


#ifndef INCLUDED_protocols_rosetta_scripts_util_hh
#define INCLUDED_protocols_rosetta_scripts_util_hh

// Unit headers

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/types.hh>
#include <basic/datacache/DataMap.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pack/task/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
// Utillity Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ headers
#include <string>
#include <set>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace rosetta_scripts {

utility::vector1< core::pack::task::operation::TaskOperationOP >
get_task_operations( utility::tag::TagCOP tag, basic::datacache::DataMap const & data );

core::pack::task::TaskFactoryOP
parse_task_operations( utility::tag::TagCOP tag, basic::datacache::DataMap const & data );

/// allows the transfer of whole taskfactories on the datamap. This way a "base" taskfactory can be created, transferred on the datamap, and
/// individual mover's specific taskoperations can be added on top
core::pack::task::TaskFactoryOP
parse_task_operations( utility::tag::TagCOP tag, basic::datacache::DataMap /*const*/ & data, core::pack::task::TaskFactoryOP & task_factory );

core::pack::task::TaskFactoryOP
parse_task_operations( std::string const & task_list, basic::datacache::DataMap const & data );

/// @brief returns a residue selector given a tag and datamap
/// @details Looks for "residue_selector" option in tag
///          If that option isn't found, returns NULL ptr
///          If that option is found, calls get_residue_selector()
core::pack::task::residue_selector::ResidueSelectorCOP
parse_residue_selector( utility::tag::TagCOP tag, basic::datacache::DataMap const & data );

/// @brief returns a residue selector given a selector's name and datamap
/// @details Looks for selector in the datamap
///          Returns a const ptr to the selector
/// @throws utility::excn::EXCN_Msg_Exception if selector is not found in datamap
core::pack::task::residue_selector::ResidueSelectorCOP
get_residue_selector( std::string const & selector_name, basic::datacache::DataMap const & data );

/// @brief Look up the score function defined in the <SCOREFXNS/>
///through the given option. Default to 'talaris2013' by default.
core::scoring::ScoreFunctionOP
parse_score_function(
	utility::tag::TagCOP tag,
	std::string const & option_name,
	basic::datacache::DataMap const & data,
	std::string const & dflt_key="talaris2013" );

/// @brief Look up the score function defined in the <SCOREFXNS/>
///through the option 'scorefxn='. Default to 'talaris2013' by default.
core::scoring::ScoreFunctionOP
parse_score_function(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data,
	std::string const & dflt_key="talaris2013" );

/// @brief Look up the name of assigned score function to the given
///option. Use this to prevent hard coding default score functions into
///protocols.
std::string
get_score_function_name(
	utility::tag::TagCOP tag,
	std::string const & option_name);

/// @brief Look up the name of assigned score function to the 'scorefxn='
///option. Use this to prevent hard coding default score functions into
///protocols.
std::string
get_score_function_name(
	utility::tag::TagCOP tag);

/// @brief convenience function to access pointers to poses that will be stored
/// in the data map at an arbitrary point during an RS protocol
/// Will look for tag in in_tag variable
core::pose::PoseOP
saved_reference_pose(
	utility::tag::TagCOP in_tag,
	basic::datacache::DataMap & data_map,
	std::string const & tag_str="reference_name" );

/// @brief variant of parse_movemap that takes in a datamap and searches it for already existing movemaps
/// Still resets movemap if MoveMap branch not found!
void
parse_movemap(
	utility::tag::TagCOP in_tag,
	core::pose::Pose const & pose,
	core::kinematics::MoveMapOP & mm,
	basic::datacache::DataMap &,
	bool const reset_movemap = true /* should we turn everything to true at start?*/);

/// @details modifies an existing movemap according to tag
/// the movemap defaults to move all bb, chi, and jumps.
/// Still resets movemap if MoveMap branch not found!
void
parse_movemap(
	utility::tag::TagCOP in_tag,
	core::pose::Pose const & pose,
	core::kinematics::MoveMapOP mm,
	bool const reset_movemap = true /* should we turn everything to true at start?*/);

/// @brief Parses in_tag, adding any MoveMaps specified in branches with names to the datamap for use after.
/// Skips any mm names that are already loaded.
/// This should allow multiple MoveMaps to be specified and loaded in a mover.
void
add_movemaps_to_datamap(
	utility::tag::TagCOP in_tag,
	core::pose::Pose const & pose,
	basic::datacache::DataMap & data,
	bool initialize_mm_as_true = false);

/// @brief Does the tag have a branch of the given name
bool
has_branch(utility::tag::TagCOP in_tag, std::string const & branch_name);

protocols::filters::FilterOP
parse_filter( std::string const & filter_name, protocols::filters::Filters_map const & d );

protocols::moves::MoverOP
parse_mover( std::string const & mover_name, protocols::moves::Movers_map const & d );

numeric::xyzVector< core::Real >
parse_xyz_vector( utility::tag::TagCOP xyz_vector_tag );


/// @brief find source residue that is nearest to res on source. If distance is greater than 2.0A, return 0. chain=0, search all chains, chain=1,2,3 etc. search only that chain
core::Size
find_nearest_res( core::pose::Pose const & source, core::pose::Pose const & target, core::Size const res, core::Size const chain = 0 );

/// @brief find nearest residue and also return the distance
void
find_nearest_res( core::pose::Pose const & source, core::pose::Pose const & target, core::Size const res, core::Size & target_res, core::Real & dist, core::Size const chain = 0 );

/// @brief finds the nearest disulife to given residue on pose by searching both up and down stream to the given postion
core::Size
find_nearest_disulfide( core::pose::Pose const & pose, core::Size const res);

/// @brief returns a vector containing all the residues with a given packer state according to the TF
utility::vector1< core::Size >
residue_packer_states( core::pose::Pose const & pose, core::pack::task::TaskFactoryCOP tf, bool const designable, bool const packable/*but not designable*/ );

/// @brief Access res_num/pdb_num tag with prefix.
/// This is to allow parsing at apply time (pdb_num) instead of at parse_my_tag.
/// This allows one to use SavePoseMovers within a protocol with the correct residue from pdb_num.
/// This function is to prevent unused variable crash.
void
parse_bogus_res_tag(utility::tag::TagCOP tag, std::string const & prefix);

} // RosettaScripts
} // protocols


#endif /*INCLUDED_protocols_RosettaScripts_util_HH*/
