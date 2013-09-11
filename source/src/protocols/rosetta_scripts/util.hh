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
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu), Rocco Moretti (rmoretti@u.washington.edu), Eva-Maria Strauch (evas01@uw.edu)

#ifndef INCLUDED_protocols_rosetta_scripts_util_hh
#define INCLUDED_protocols_rosetta_scripts_util_hh

// Unit headers

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/types.hh>
#include <protocols/moves/DataMap.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
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
get_task_operations( utility::tag::TagPtr const tag, protocols::moves::DataMap const & data );

core::pack::task::TaskFactoryOP
parse_task_operations( utility::tag::TagPtr const tag, protocols::moves::DataMap const & data );

/// allows the transfer of whole taskfactories on the datamap. This way a "base" taskfactory can be created, transferred on the datamap, and
/// individual mover's specific taskoperations can be added on top
core::pack::task::TaskFactoryOP
parse_task_operations( utility::tag::TagPtr const tag, protocols::moves::DataMap /*const*/ & data, core::pack::task::TaskFactoryOP & task_factory );

core::pack::task::TaskFactoryOP
parse_task_operations( std::string const task_list, protocols::moves::DataMap const & data );

///@brief Look up the score function defined in the <SCOREFXNS/>
///through the given option. Default to 'talaris2013' by default.
core::scoring::ScoreFunctionOP
parse_score_function(
	utility::tag::TagPtr const tag,
	std::string const & option_name,
	protocols::moves::DataMap const & data,
	std::string const dflt_key="talaris2013" );

///@brief Look up the score function defined in the <SCOREFXNS/>
///through the option 'scorefxn='. Default to 'talaris2013' by default.
core::scoring::ScoreFunctionOP
parse_score_function(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap const & data,
	std::string const dflt_key="talaris2013" );

///@brief Look up the name of assigned score function to the given
///option. Use this to prevent hard coding default score functions into
///protocols.
std::string
get_score_function_name(
	utility::tag::TagPtr const tag,
	std::string const & option_name);

///@brief Look up the name of assigned score function to the 'scorefxn='
///option. Use this to prevent hard coding default score functions into
///protocols.
std::string
get_score_function_name(
	utility::tag::TagPtr const tag);

/// @brief convenience function to access pointers to poses that will be stored
/// in the data map at an arbitrary point during an RS protocol
core::pose::PoseOP
saved_reference_pose( utility::tag::TagPtr const in_tag, protocols::moves::DataMap & data_map );

void
parse_movemap( utility::tag::TagPtr const in_tag, core::pose::Pose const & pose, core::kinematics::MoveMapOP & mm, protocols::moves::DataMap &, bool const reset_movemap = true /* should we turn everything to true at start?*/ );

void
parse_movemap( utility::tag::TagPtr const in_tag, core::pose::Pose const & pose, core::kinematics::MoveMapOP mm );

protocols::filters::FilterOP
parse_filter( std::string const filter_name, protocols::filters::Filters_map const & d );

protocols::moves::MoverOP
parse_mover( std::string const mover_name, protocols::moves::Movers_map const & d );

numeric::xyzVector< core::Real >
parse_xyz_vector( utility::tag::TagPtr const xyz_vector_tag );


/// @brief find source residue that is nearest to res on source. If distance is greater than 2.0A, return 0. chain=0, search all chains, chain=1,2,3 etc. search only that chain
core::Size
find_nearest_res( core::pose::Pose const & source, core::pose::Pose const & target, core::Size const res, core::Size const chain = 0 );

/// @brief returns a vector containing all the residues with a given packer state according to the TF
utility::vector1< core::Size >
residue_packer_states( core::pose::Pose const & pose, core::pack::task::TaskFactoryCOP tf, bool const designable, bool const packable/*but not designable*/ );

} // RosettaScripts
} // protocols


#endif /*INCLUDED_protocols_RosettaScripts_util_HH*/
