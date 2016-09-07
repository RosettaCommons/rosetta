// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Neil King (neilking@uw.edu)

#ifndef INCLUDED_devel_matdes_StoreQuasiSymmetricTaskMover_fwd_hh
#define INCLUDED_devel_matdes_StoreQuasiSymmetricTaskMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace matdes {

class StoreQuasiSymmetricTaskMover;

//typedef utility::pointer::owning_ptr< StoreQuasiSymmetricTaskMover > StoreQuasiSymmetricTaskMoverOP;
//typedef utility::pointer::owning_ptr< StoreQuasiSymmetricTaskMover const > StoreQuasiSymmetricTaskMoverCOP;
typedef utility::pointer::shared_ptr< StoreQuasiSymmetricTaskMover > StoreQuasiSymmetricTaskMoverOP;
typedef utility::pointer::shared_ptr< StoreQuasiSymmetricTaskMover const > StoreQuasiSymmetricTaskMoverCOP;

} // matdes
} // devel

#endif
