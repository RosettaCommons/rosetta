// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/RotamerSetOperations/SpecialRotamerRotSetOps.fwd.hh
/// @brief  combine two rotamer sets
/// @author Summer Thyme, sthyme@u.washington.edu

#ifndef INCLUDED_protocols_toolbox_rotamer_set_operations_SpecialRotamerRotSetOps_fwd_hh
#define INCLUDED_protocols_toolbox_rotamer_set_operations_SpecialRotamerRotSetOps_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {
namespace rotamer_set_operations {

class SpecialRotamerRSO;

typedef utility::pointer::shared_ptr< SpecialRotamerRSO > SpecialRotamerRSOOP;
typedef utility::pointer::shared_ptr< SpecialRotamerRSO const > SpecialRotamerRSOCOP;

} //namespace rotamer_set_operations
} //namespace toolbox
} //namespace protocols

#endif
