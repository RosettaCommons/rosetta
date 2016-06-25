// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/rna/leontis_westhof_util.hh
/// @brief  Implementation of Leontis/Westhof nucleic acid base-pair classification
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_rna_leontis_westhof_util_HH
#define INCLUDED_core_pose_rna_leontis_westhof_util_HH

#include <core/chemical/rna/util.hh>

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core::chemical::rna;

namespace core {
namespace pose {
namespace rna {

LW_BaseDoubletOrientation
get_LW_orientation( BaseEdge const & edge1, BaseEdge const & edge2, BaseDoubletOrientation const & orientation );

BaseDoubletOrientation
get_base_doublet_orientation_from_LW( BaseEdge const & edge1, BaseEdge const & edge2, LW_BaseDoubletOrientation const & lw_orientation );

typedef std::map< std::pair< BaseEdge, BaseEdge >, std::map< BaseDoubletOrientation, LW_BaseDoubletOrientation > > LW_Table;

LW_Table const &
get_leontis_westhof_table();

} //rna
} //pose
} //core

#endif
