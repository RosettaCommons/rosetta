// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/MPI_WorkUnitManager.fwd.hh
/// @brief  MPI_WorkUnitManager base classes
/// @author Mike Tyka

#ifndef INCLUDED_protocols_wum_MPI_WorkUnitManager_fwd_hh
#define INCLUDED_protocols_wum_MPI_WorkUnitManager_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace wum {

// Forward
class MPI_WorkUnitManager;

typedef utility::pointer::shared_ptr< MPI_WorkUnitManager > MPI_WorkUnitManagerOP;
typedef utility::pointer::shared_ptr< MPI_WorkUnitManager const > MPI_WorkUnitManagerCOP;

} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_MPI_WorkUnitManager_FWD_HH

