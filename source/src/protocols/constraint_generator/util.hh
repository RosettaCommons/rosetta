// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/constraint_generator/util.hh
/// @brief Utility functions for constraint generators
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_constraint_generator_util_hh
#define INCLUDED_protocols_constraint_generator_util_hh

#include <core/id/SequenceMapping.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace constraint_generator {

/// @brief gets native pose from command line option, if one is specified
core::pose::PoseOP
get_native_pose();

/// @brief helper function to wrap a function with a scalar weighted function
core::scoring::func::FuncOP
scalar_weighted( core::scoring::func::FuncOP func, core::Real const weight );

/// @brief creates a function from a text definition
core::scoring::func::FuncOP
create_func( std::string const & func_def );

/// @brief creates a sequencemapping from pose1 to pose2
core::id::SequenceMapping
generate_seqmap_from_poses(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2 );

/// @brief find number of residues -- if pose is symetric,
/// returns nubmer of symmetry-independent residues
core::Size
compute_nres_in_asymmetric_unit( core::pose::Pose const & pose );

} //protocols
} //constraint_generator


#endif	//protocols/constraint_generator_util_hh

