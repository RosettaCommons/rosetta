// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Align a random jump to template
/// @detailed
/// @author Yifan Song

#ifndef INCLUDED_protocols_moves_util_hh
#define INCLUDED_protocols_moves_util_hh

#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Edge.hh>

#include <core/types.hh>

#include <list>

namespace protocols {
namespace comparative_modeling {
namespace hybridize {
			
using namespace core;
using namespace kinematics;

bool discontinued_upper(core::pose::Pose const & pose, Size const seqpos);

bool discontinued_lower(core::pose::Pose const & pose, Size const seqpos);

std::list < Size > downstream_residues_from_jump(core::pose::Pose const & pose, Size const jump_number);

} // hybridize 
} // comparative_modeling 
} // protocols

#endif

