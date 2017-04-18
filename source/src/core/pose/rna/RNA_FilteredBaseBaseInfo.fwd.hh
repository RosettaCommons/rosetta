// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/rna/RNA_FilteredBaseBasePotential.fwd.hh
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das

#ifndef INCLUDED_core_pose_rna_RNA_FilteredBaseBaseInfo_fwd_hh
#define INCLUDED_core_pose_rna_RNA_FilteredBaseBaseInfo_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pose {
namespace rna {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of RNA centroid information inside the pose.
class RNA_FilteredBaseBaseInfo;
typedef utility::pointer::shared_ptr< RNA_FilteredBaseBaseInfo > RNA_FilteredBaseBaseInfoOP;

} //rna
} //scoring
} //core

#endif
