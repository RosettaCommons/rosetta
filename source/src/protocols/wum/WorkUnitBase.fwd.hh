// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/WorkUnitBase.fwd.hh
/// @brief  WorkUnitBase base classes
/// @author Mike Tyka

#ifndef INCLUDED_protocols_wum_WorkUnitBase_fwd_hh
#define INCLUDED_protocols_wum_WorkUnitBase_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace wum {

class WorkUnitBase;
typedef utility::pointer::shared_ptr< WorkUnitBase > WorkUnitBaseOP;
typedef utility::pointer::shared_ptr< WorkUnitBase const > WorkUnitBaseCOP;

class WorkUnit_Wait;
typedef utility::pointer::shared_ptr< WorkUnit_Wait > WorkUnit_WaitOP;
typedef utility::pointer::shared_ptr< WorkUnit_Wait const > WorkUnit_WaitCOP;

class WorkUnit_SilentStructStore;
typedef utility::pointer::shared_ptr< WorkUnit_SilentStructStore > WorkUnit_SilentStructStoreOP;
typedef utility::pointer::shared_ptr< WorkUnit_SilentStructStore const > WorkUnit_SilentStructStoreCOP;

class WorkUnit_MoverWrapper;
typedef utility::pointer::shared_ptr< WorkUnit_MoverWrapper > WorkUnit_MoverWrapperOP;
typedef utility::pointer::shared_ptr< WorkUnit_MoverWrapper const > WorkUnit_MoverWrapperCOP;

} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_WorkUnitBase_FWD_HH

