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
#include <core/pose/Pose.fwd.hh>


namespace protocols {
namespace toolbox {
namespace match_enzdes_util{

void
replace_residue_keeping_all_atom_positions(
	core::pose::Pose & pose,
	core::conformation::Residue new_res,
	core::Size res_pos
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

}
} // enzdes
} //protocols


#endif
