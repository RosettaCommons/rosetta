// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file IO-functionality for enzyme Constraints
/// @brief
/// @author Florian Richter, floric@u.washington.edu

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_util_functions_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_util_functions_hh



#include <core/types.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/constraints/AmbiguousConstraint.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>

#include <list>

namespace protocols {
namespace toolbox {
namespace match_enzdes_util{

void
replace_residue_keeping_all_atom_positions(
	core::pose::Pose & pose,
	core::conformation::Residue new_res,
	core::Size res_pos
);

/// @brief fowarding function for function below
core::scoring::constraints::AmbiguousConstraintCOP
constrain_pose_res_to_invrots(
	std::list< core::conformation::ResidueCOP> const & invrots,
	utility::vector1< core::Size > const & seqpos,
	core::pose::Pose const & pose,
	core::scoring::constraints::FuncOP constraint_func = NULL
);

/// @brief constraints each invrot to the
/// backbone of each seqpos and throws all
/// those constraints into one ambiguous
/// constraint.
core::scoring::constraints::AmbiguousConstraintCOP
constrain_pose_res_to_invrots(
	std::list< core::conformation::ResidueCOP> const & invrots,
	utility::vector1< core::Size > const & seqpos,
	core::pose::Pose const & pose,
	core::id::AtomID const & fixed_pt,
	core::scoring::constraints::FuncOP constraint_func = NULL
);


/// @brief convenience function that returns a residue
/// of the desired cst interaction
/// in case there are no constraints in the pose,
/// returns null pointer
core::conformation::ResidueCOP
cst_residue_in_pose(
	core::pose::Pose const & pose,
	core::Size geomcst,
	core::Size geomcst_template_res
);

std::string
assemble_remark_line(
	std::string chainA,
	std::string resA,
	int seqposA,
	std::string chainB,
	std::string resB,
	int seqposB,
	core::Size cst_block,
	core::Size ex_geom_id = 1
);


bool
split_up_remark_line(
	std::string line,
	std::string & chainA,
	std::string & resA,
	int & seqposA,
	std::string & chainB,
	std::string & resB,
	int & seqposB,
	core::Size & cst_block,
	core::Size & ex_geom_id
);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// @brief finds the first non-ligand residue in the pose  (should be the N-terminus)
core::Size get_first_protein_residue( core::pose::Pose const & pose );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// @brief finds the first non-ligand residue in the pose  (should be the N-terminus)
core::Size get_last_protein_residue( core::pose::Pose const & pose );


}  // match_enzdes_util
} // toolbox
} //protocols


#endif
