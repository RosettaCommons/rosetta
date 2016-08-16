// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MessageListener.fwd.hh
///
/// @brief
/// @author Tim Jacobs


#ifndef INCLUDED_basic_mpi_MessageListener_FWD_HH
#define INCLUDED_basic_mpi_MessageListener_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace basic {
namespace mpi {

/// @brief MessageListener id tags
///note: please add string conversion to util.cc/hh
enum listener_tags {
	DATABASE_PROTOCOL_AND_BATCH_ID_TAG = 42
};


class MessageListener;
typedef utility::pointer::shared_ptr< MessageListener > MessageListenerOP;
typedef utility::pointer::shared_ptr< MessageListener const > MessageListenerCOP;

} //namespace
} //namespace
#endif


