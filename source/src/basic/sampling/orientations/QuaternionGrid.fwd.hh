// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/sampling/orientations/QuaternionGrid.fwd.hh
/// @brief  basic::sampling::orientations::QuaternionGrid forward declaration
/// @author Will Sheffler <willsheffler@gmail.com>

#ifndef INCLUDED_basic_sampling_orientations_QuaternionGrid_fwd_hh
#define INCLUDED_basic_sampling_orientations_QuaternionGrid_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace basic {
namespace sampling {
namespace orientations {

class QuaternionGrid;
typedef utility::pointer::shared_ptr< QuaternionGrid const > QuaternionGridCOP;
typedef utility::pointer::weak_ptr< QuaternionGrid const > QuaternionGridCAP;

class QuaternionGridManager;
typedef utility::pointer::weak_ptr< QuaternionGridManager const > QuaternionGridManagerCAP;

}
}
}


#endif /*INCLUDED_basic_sampling_orientations_QuaternionGrid_FWD_HH*/
