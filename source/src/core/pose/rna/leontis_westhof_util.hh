// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/rna/leontis_westhof_util.hh
/// @brief  Implementation of Leontis/Westhof nucleic acid base-pair classification
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_rna_leontis_westhof_util_HH
#define INCLUDED_core_pose_rna_leontis_westhof_util_HH

#include <core/chemical/rna/util.hh>

namespace core {
namespace pose {
namespace rna {

core::chemical::rna::LW_BaseDoubletOrientation
get_LW_orientation(
	core::chemical::rna::BaseEdge const & edge1,
	core::chemical::rna::BaseEdge const & edge2,
	core::chemical::rna::BaseDoubletOrientation const & orientation
);

core::chemical::rna::BaseDoubletOrientation
get_base_doublet_orientation_from_LW(
	core::chemical::rna::BaseEdge const & edge1,
	core::chemical::rna::BaseEdge const & edge2,
	core::chemical::rna::LW_BaseDoubletOrientation const & lw_orientation
);

typedef std::map< std::pair< core::chemical::rna::BaseEdge, core::chemical::rna::BaseEdge >, std::map< core::chemical::rna::BaseDoubletOrientation, core::chemical::rna::LW_BaseDoubletOrientation > > LW_Table;

LW_Table const &
get_leontis_westhof_table();

} //rna
} //pose
} //core

#endif
