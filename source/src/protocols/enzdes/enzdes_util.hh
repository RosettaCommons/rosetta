// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/enzdes/enzdes_util.hh
/// @brief a bunch of utility functions used in enzdes
/// @author Florian Richter, floric@u.washington.edu


#ifndef INCLUDED_protocols_enzdes_enzdes_util_hh
#define INCLUDED_protocols_enzdes_enzdes_util_hh


// Unit headers


// Package headers


// Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
//for parse_task_operations

#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.fwd.hh>
// Utility Headers

#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>


//Utility Headers

// C++ Headers

namespace protocols {
namespace enzdes {
namespace enzutil {

/// @brief queries whether a certain position is considered catalytic,
/// i.e. if the residue is constrained according to the enzdes constraint
/// file
bool
is_catalytic_seqpos(
	core::pose::Pose const & pose,
	core::Size const seqpos
);

utility::vector1< core::Size >
catalytic_res( core::pose::Pose const & pose);

/// @brief sums up coordinate, atom_pair, angle, and dihedral constraint
/// scores for the residue in question. in case which_res is -1, the scores
/// for the whole pose will be taken
core::Real
sum_constraint_scoreterms(
	core::pose::Pose const & pose,
	int which_res
);


void
read_pose_from_file(
	core::pose::Pose & pose,
	std::string const & filename
);

void
make_continuous_true_regions_in_bool_vector(
	utility::vector1< bool > & the_vector,
	core::Size const min_number_continuous_trues
);

core::pack::task::PackerTaskOP
recreate_task(
	core::pose::Pose const & pose,
	core::pack::task::PackerTask const & orig_task
);

/// @brief convenience function to get enzdes cst io out of
/// cst cache. note: may return empty pointer in case
/// cst cache wasn't set yet
toolbox::match_enzdes_util::EnzConstraintIOCOP
get_enzcst_io( core::pose::Pose const & pose );

void
remove_remark_header_for_geomcst(
	core::pose::Pose & pose,
	core::Size geomcst
);

/// @brief function to write proper remark headers
/// from whatever is found in the pose cstcache
void
create_remark_headers_from_cstcache(
	core::pose::Pose & pose
);

std::string
get_pdb_code_from_pose_tag( core::pose::Pose const & pose );

/// @TODO delete and use cctype functions
bool
is_digit( char * cha );

/// @TODO delete and use cctype functions
bool
is_uppercase_letter( char * cha);

/// @TODO delete and use cctype functions
bool
is_lowercase_letter( char * cha );

void
disable_constraint_scoreterms(core::scoring::ScoreFunctionOP scorefxn);
void
enable_constraint_scoreterms(core::scoring::ScoreFunctionOP scorefxn);
bool
is_scofx_cstfied(core::scoring::ScoreFunctionCOP scorefxn);
void
scorefxn_update_from_options(core::scoring::ScoreFunctionOP scorefxn);

void
remove_all_enzdes_constraints( core::pose::Pose & pose );

void
get_resnum_from_cstid_list( std::string const& cstidlist, core::pose::Pose const& pose, utility::vector1<core::Size>& resnums );

core::Size
get_resnum_from_cstid( std::string const& cstid, core::pose::Pose const & pose);


}//enzutil
} //enzdes
} //protocols


#endif
