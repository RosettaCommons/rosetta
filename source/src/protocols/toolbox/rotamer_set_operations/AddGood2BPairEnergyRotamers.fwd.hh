// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/RotamerSetOperations/AddGood2BPairEnergyRotamers.fwd.hh
/// @brief
/// @author Florian Richter

#ifndef INCLUDED_protocols_toolbox_rotamer_set_operations_AddGood2BPairEnergyRotamers_fwd_hh
#define INCLUDED_protocols_toolbox_rotamer_set_operations_AddGood2BPairEnergyRotamers_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {
namespace rotamer_set_operations {

class AddGood2BPairEnergyRotamers;

typedef utility::pointer::shared_ptr< AddGood2BPairEnergyRotamers > AddGood2BPairEnergyRotamersOP;
typedef utility::pointer::shared_ptr< AddGood2BPairEnergyRotamers const > AddGood2BPairEnergyRotamersCOP;

} //namespace rotamer_set_operations
} //namespace toolbox
} //namespace protocols

#endif
