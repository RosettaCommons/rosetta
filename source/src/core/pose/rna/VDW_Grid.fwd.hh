// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/rna/VDW_Grid.fwd.hh
/// @brief
/// @details
/// @author Caleb Geniesse, geniesse@stanford.edu


#ifndef INCLUDED_core_pose_rna_VDW_Grid_FWD_HH
#define INCLUDED_core_pose_rna_VDW_Grid_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pose {
namespace rna {

struct Atom_Bin;

class VDW_Grid;
typedef utility::pointer::shared_ptr< VDW_Grid > VDW_GridOP;
typedef utility::pointer::shared_ptr< VDW_Grid const > VDW_GridCOP;

} //rna
} //pose
} //core

#endif
